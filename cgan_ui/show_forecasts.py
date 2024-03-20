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
from datetime import timedelta
import cfgrib
import xarray as xr


# Required information about each variable we want to look at
var_info = {
    "sp": {
        "name": "Surface pressure",
        "units": "hPa",
        "normalisation": 100,  # Convert from Pa to hPa
        "offset": 0,
        "accumulated": False,
    },
    "msl": {
        "name": "Pressure at mean sea level",
        "units": "hPa",
        "normalisation": 100,  # Convert from Pa to hPa
        "offset": 0,
        "accumulated": False,
    },
    "t2m": {
        "name": "Two metre temperature",
        "units": "deg. C",
        "normalisation": 1,
        "offset": -273.15,  # Convert from Kelvin to deg. C
        "accumulated": False,
    },
    "wind": {
        "name": "Wind speed",
        "units": "m/s",
        "normalisation": 1,
        "offset": 0,
        "accumulated": False,
    },
    "tp": {
        "name": "Total precipitation",
        "units": "mm/day",
        "normalisation": 0.001,  # Convert from m to mm/day
        "offset": 0,
        "accumulated": True,
    },
    "ro": {
        "name": "Surface runoff water",
        "units": "m",
        "normalisation": 1,
        "offset": 0,
        "accumulated": False,
    },
}


# All forecasts are initialised at 00:00 UTC
# Lead times are 30, 33, 36, 39, 42, 45, 48, 51, 54 hours.
start_hour = 30
end_hour = 54


# Are all dimensions present in a data set
#    dataset    - An xarray Dataset or xarray DataArray
#    dimensions - A string representing a dimension name or a list of strings
#                 representing multiple dimension names.
# Returns
#    Either True of False
def all_dimensions_present(dataset, dimensions):

    # Is dataset an xarray Dataset or an xarray DataArray
    if (type(dataset) is xr.core.dataarray.Dataset) or (
        type(dataset) is xr.core.dataarray.DataArray
    ):

        # If no dimensions are specified
        if dimensions is None:
            return True

        # If dimensions is a single string
        elif type(dimensions) is str:
            # Is this dimension one of the dimensions?
            return dimensions in list(dataset.dims)

        # Is coordinates a list
        elif type(dimensions) is list:
            # Before we start assume all dimensions are present
            all_present = True
            # For each list element
            for dim in dimensions:
                # Is this coord and all previous coordinates present?
                all_present = all_present and (dim in list(dataset.dims))
            return all_present

        # dimensions is not a list or a single string
        else:
            return False

    # dataset is not an xarray Dataset or an xarray DataArray
    else:
        return False


# Find the first xarray DataArray corresponding to a key somewhere in a list of xarray Datasets
#    dataset_list - list of xarray Datasets or an xarray Dataset or an xarray DataArray
#    key          - string corresponding to a variable in one of the datasets
#    dimensions   - Optional string or list of strings corresponding to dimensions that must be present
# Returns the corresponding DataArray if it exists or None if  it doesn't.
def DataArray_from_Dataset_list(dataset_list, key, dimensions=None):

    # Is dataset_list an xarray DataArray?
    if type(dataset_list) is xr.core.dataarray.DataArray:
        # Does the DataArray have the same name as key and all necessary dimensions?
        if (dataset_list.name == key) and all_dimensions_present(
            dataset_list, dimensions
        ):
            # This is the correct DataArray
            return dataset_list

    # Is dataset_list an xarray Dataset?
    elif type(dataset_list) is xr.core.dataarray.Dataset:
        # The keys in the Dataset
        keys = list(dataset_list.keys())
        # If key is in the Dataset and all dimensions are present
        if key in keys and all_dimensions_present(dataset_list[key], dimensions):
            # Return the DataArray corresponding to key
            return dataset_list[key]

    # Is dataset_list a list?
    elif type(dataset_list) is list:
        # For each element in the list
        for i in range(len(dataset_list)):
            # Is dataset_list element an xarray Dataset?
            if type(dataset_list[i]) is xr.core.dataarray.Dataset:
                # The keys in the Dataset
                keys = list(dataset_list[i].keys())
                # If key is in the Dataset and all dimensions are present
                if key in keys and all_dimensions_present(dataset_list[i], dimensions):
                    # Return the DataArray corresponding to key
                    return dataset_list[i][key]


# Load a 24 hour mean forecast at a lead time of 30 to 54 hours
# To be clear about the dates use print_forecast_info(forecast_init_date)
# Arguments
#    key                  - The variable to load. Options are returned by get_possible_variables()
#    forecast_init_date   - A datetime.datetime corresponding to when the forecast was initialised.
#    data_dir             - Directory where the data is stored. Same filename format as used at ECMWF.
#    status_updates=True  - Show what files are currently loading (True or False)
# Returns
#    An xarray DataArray containing the 24 hour mean quantity.
def load_forecast(key, forecast_init_date, data_dir, status_updates=True):

    # Shorthand for clarity
    d = forecast_init_date

    if var_info[key]["accumulated"]:  # Accumulated variables
        forecast_hours = [start_hour, end_hour]
        if start_hour == end_hour:
            print(
                "Error in load_forecast, start_hour == end_hour for accumulated forecast."
            )

    else:  # Variables are not accumulated
        forecast_hours = np.arange(start_hour, end_hour + 1, 3)

    # Just need the start and end lead times for accumulated variables
    for lead_hour in forecast_hours:
        if status_updates:
            print(f"Loading {key} with lead time {lead_hour}h.")

        # Name of the file we will read
        file_name = f"{data_dir}/{d.year}{d.month:02}{d.day:02}000000-{lead_hour}h-enfo-ef.grib2"

        # Open a grib2 file for reading
        # xarray.open_dataset(file_name, engine="cfgrib") doesn't work!
        ds = cfgrib.open_datasets(file_name)

        # Get the DataArray corresponding to the key. 'number' ensures that we pick the ensemble forecast.
        da = DataArray_from_Dataset_list(ds, key, dimensions="number")

        # Get the DataArray corresponding to the key. 'number' ensures that we pick the ensemble forecast.
        if key == "wind":
            da_u10 = DataArray_from_Dataset_list(ds, "u10", dimensions="number")
            da_v10 = DataArray_from_Dataset_list(ds, "u10", dimensions="number")
            # We want the average wind speed, not the average air velocity.
            data_3h = np.sqrt(
                da_u10[:, 261:416, 796:938] ** 2 + da_v10[:, 261:416, 796:938] ** 2
            )
        else:
            da = DataArray_from_Dataset_list(ds, key, dimensions="number")
            data_3h = da[:, 261:416, 796:938]

        # Average over the forecast period
        if lead_hour == start_hour:
            data = data_3h / 2  # Divide by two for the trapezium rule
            if var_info[key]["accumulated"]:
                data = -data  # Negative because we are taking it away
            data_norm = 0.5

        elif lead_hour == end_hour:
            data += data_3h / 2  # Divide by two for the trapezium rule
            data_norm += 0.5

        else:
            data += data_3h
            data_norm += 1

    # Normalise the data to get the 24 hour mean with the correct units
    data = var_info[key]["offset"] + data / (data_norm * var_info[key]["normalisation"])

    # Record the name and units
    data.attrs["name"] = var_info[key]["name"]
    data.attrs["units"] = var_info[key]["units"]

    return data


# Plot the ensemble mean and ensemble standard deviation of the forecast data
def plot_forecast(data):

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
        data["longitude"],
        data["latitude"],
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
        data["longitude"],
        data["latitude"],
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


# Some info to be clear with dates and times
# Arguments
#    forecast_init_date - A datetime.datetime corresponding to when the forecast was initialised.
def print_forecast_info(forecast_init_date):
    start_date = forecast_init_date + timedelta(hours=start_hour)
    end_date = forecast_init_date + timedelta(hours=end_hour)
    print(f"Forecast average: {start_date} - {end_date}")
    print(f"Forecast initialisation: {forecast_init_date.date()} 00:00:00")
    print()


# Lists the variables available for plotting
# Arguments
#    print_variables=True - print a list of the variables (True or False)
# Returns
#    A list of the variable names that can be used for plotting.
def get_possible_variables(print_variables=True):

    keys = list(var_info.keys())
    if print_variables:
        print("Available variables to plot are the following:")
        for i in range(len(keys)):
            print(
                f"{keys[i]:<5} - {var_info[keys[i]]['name']} ({var_info[keys[i]]['units']})"
            )
        print()

    return list(var_info.keys())
