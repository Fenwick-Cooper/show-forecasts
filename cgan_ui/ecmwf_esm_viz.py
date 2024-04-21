# Load and plot the ECMWF open forecast data

# To do:
#    Option to specify forecast lead times
#       In particular specify 7 and 14 day forecasts
#    Arbitrary region specification
#    Plot relative humidity
#       Requires interpolating in the vertical to the correct surface pressure
#    Plot ensemble members
#    Runoff colour scale is not the best
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from datetime import datetime
import xarray as xr
from cgan_ui.constants import DATA_PARAMS, LEAD_START_HOUR, LEAD_END_HOUR


# Load a 24 hour mean forecast at a lead time of 30 to 54 hours
# To be clear about the dates use print_forecast_info(forecast_init_date)
# Arguments
#    key                  - The variable to load. Options are returned by get_possible_variables()
#    forecast_init_date   - A datetime.datetime corresponding to when the forecast was initialised.
#    data_dir             - Directory where the data is stored. Same filename format as used at ECMWF.
#    status_updates=True  - Show what files are currently loading (True or False)
#    file_ext=nc          - forecast dataset file extension. defaults to nc (NetCDF)
# Returns
#    An xarray DataArray containing the 24 hour mean quantity.
def load_forecast(
    key: str,
    forecast_init_date: datetime.date,
    data_dir: str,
    status_updates: bool | None = True,
    file_ext: str | None = "nc",
):

    # Shorthand for clarity
    d = forecast_init_date

    if DATA_PARAMS[key]["accumulated"]:  # Accumulated variables
        forecast_hours = [LEAD_START_HOUR, LEAD_END_HOUR]
        if LEAD_START_HOUR == LEAD_END_HOUR:
            print(
                "Error in load_forecast, start_hour == end_hour for accumulated forecast."
            )

    else:  # Variables are not accumulated
        forecast_hours = np.arange(LEAD_START_HOUR, LEAD_END_HOUR + 1, 3)

    # Just need the start and end lead times for accumulated variables
    for lead_hour in forecast_hours:
        # Name of the file we will read
        file_name = f"{data_dir}/{d.year}{d.month:02}{d.day:02}000000-{lead_hour}h-enfo-ef.{file_ext}"

        if status_updates:
            print(
                f"Loading {key} with lead time {lead_hour}h from {file_name.split('/')[-1]}"
            )

        # Open a NetCDF file for reading
        ds = xr.open_dataset(file_name)

        # Get the DataArray corresponding to the key. 'number' ensures that we pick the ensemble forecast.
        if key == "wind":
            # We want the average wind speed, not the average air velocity.
            data_3h = np.sqrt(ds["u10"] ** 2 + ds["v10"] ** 2)
        else:
            data_3h = ds[key]

        # Average over the forecast period
        if lead_hour == LEAD_START_HOUR:
            data = data_3h / 2  # Divide by two for the trapezium rule
            if DATA_PARAMS[key]["accumulated"]:
                data = -data  # Negative because we are taking it away
            data_norm = 0.5

        elif lead_hour == LEAD_END_HOUR:
            data += data_3h / 2  # Divide by two for the trapezium rule
            data_norm += 0.5

        else:
            data += data_3h
            data_norm += 1

    # Normalise the data to get the 24 hour mean with the correct units
    data = DATA_PARAMS[key]["offset"] + data / (
        data_norm * DATA_PARAMS[key]["normalisation"]
    )

    # Record the name and units
    data.attrs["name"] = DATA_PARAMS[key]["name"]
    data.attrs["units"] = DATA_PARAMS[key]["units"]

    return data


# Plot the ensemble mean and ensemble standard deviation of the forecast data
def plot_forecast(
    data: xr.DataArray, lon_dim: str | None = "lon", lat_dim: str | None = "lat"
):

    # Define the figure and each axis for the rows and columns
    fig, axs = plt.subplots(
        nrows=1, ncols=2, subplot_kw={"projection": ccrs.PlateCarree()}, figsize=(10, 5)
    )

    # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs = axs.flatten()

    ax = axs[0]  # First plot (left)
    ax.add_feature(
        cfeature.COASTLINE, linewidth=1
    )  # Draw some features to see where we are
    ax.add_feature(cfeature.BORDERS, linewidth=1)
    ax.add_feature(
        cfeature.LAKES,
        linewidth=1,
        linestyle="-",
        edgecolor="dimgrey",
        facecolor="none",
    )
    # Actually make the plot
    c = ax.pcolormesh(
        data[lon_dim],
        data[lat_dim],
        np.mean(data, axis=0),
        cmap="YlGnBu",
        transform=ccrs.PlateCarree(),
    )
    cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
    cb.set_label(data.attrs["units"])  # Label the colorbar
    ax.set_title(f"Ensemble mean", size=14)  # This plot's title

    ax = axs[1]  # Second plot (right)
    ax.add_feature(
        cfeature.COASTLINE, linewidth=1
    )  # Draw some features to see where we are
    ax.add_feature(cfeature.BORDERS, linewidth=1)
    ax.add_feature(
        cfeature.LAKES,
        linewidth=1,
        linestyle="-",
        edgecolor="dimgrey",
        facecolor="none",
    )
    # Actually make the plot
    c = ax.pcolormesh(
        data[lon_dim],
        data[lat_dim],
        np.std(data, axis=0, ddof=1),
        cmap="YlGnBu",
        transform=ccrs.PlateCarree(),
    )
    cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
    cb.set_label(data.attrs["units"])  # Label the colorbar
    ax.set_title(f"Ensemble standard deviation", size=14)  # This plot's title

    fig.suptitle(data.attrs["name"])  # Overall title
    plt.tight_layout()  # Looks nicer
    plt.show()  # Finally draw the plot
