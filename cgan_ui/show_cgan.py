# Load and plot cGAN forecast data.

# To do:
#   Change "mm h**-1" to "mm/h" in the data.
#   Store the model name in the data.
#   Combine with the show_forecasts.py script used for plotting IFS open data.

import numpy as np
from pathlib import Path
from os import getenv
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import colors  # For consistency with Harris et. al 2022
from datetime import datetime, timedelta
import xarray as xr
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
from matplotlib import colors  # For consistency with Harris et. al 2022
import xarray as xr
from cgan_ui.constants import COUNTRY_NAMES, COLOR_SCHEMES
from cgan_ui.data_utils import (
    get_contour_levels,
    get_plot_normalisation,
    get_shape_boundary,
    get_region_extent,
    datetime64_to_datetime,
)


# Load a 24 hour mean forecast at a lead time of 30 to 54 hours
# To be clear about the dates use print_forecast_info(forecast_init_date)
# Arguments
#    forecast_init_date   - A datetime.datetime corresponding to when the forecast was initialised.
#    data_dir             - Directory where the data is stored.
# Returns
#    An xarray DataSet containing the cGAN rainfall forecasts.
def load_GAN_forecast(
    forecast_init_date: datetime, data_dir: Path, mask_region: str
) -> xr.Dataset:
    d = forecast_init_date  # Shorthand
    file_path = (
        data_dir
        / str(d.year)
        / str(d.month).rjust(2, "0")
        / f"{mask_region.lower().replace(' ','_')}-cgan_forecast-{d.year}{d.month:02}{d.day:02}.nc"
    )
    data = xr.open_dataset(file_path)
    return data


# Plot the ensemble mean and ensemble standard deviation of the cGAN forecast data
# at each valid time.
# Arguments
#   data              - An xarray DataSet containing the cGAN rainfall forecasts.
#   style=None        - Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
#   plot_units='mm/h' - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   region='ICPAC'    - can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                       'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
def plot_GAN_forecast(
    data: xr.Dataset,
    lon_dim: str | None = "longitude",
    lat_dim: str | None = "latitude",
    style: str | None = COLOR_SCHEMES[0],
    plot_units: str | None = "mm/h",
    region: str | None = COUNTRY_NAMES[0],
):

    # Get the units to use for plotting
    plot_norm = get_plot_normalisation(plot_units)

    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15 * plot_norm)

    # Use a style other than the default
    if style is not None:
        plot_levels, plot_colours = get_contour_levels(style)

    # Load the border shapefile
    reader = get_shape_boundary(shape_name=region)
    shape_feature = ShapelyFeature(
        reader.geometries(), ccrs.PlateCarree(), facecolor="none"
    )

    # Get the extent of the region that we are looking at
    if region != COUNTRY_NAMES[0]:
        region_extent = get_region_extent(region, border_size=0.5)

    # There are plots for each valid time
    for valid_time_idx in range(len(data["valid_time"])):

        # Convert the forecast valid time to a datetime.datetime format
        valid_time = datetime64_to_datetime(
            data["fcst_valid_time"][0, valid_time_idx].values
        )

        # Define the figure and each axis for the rows and columns
        fig, axs = plt.subplots(
            nrows=1,
            ncols=2,
            subplot_kw={"projection": ccrs.PlateCarree()},
            figsize=(10, 5),
        )

        # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
        axs = axs.flatten()

        ax = axs[0]  # First plot (left)
        ax.add_feature(
            cfeature.COASTLINE, linewidth=1
        )  # Draw some features to see where we are
        ax.add_feature(shape_feature)  # The borders
        ax.add_feature(
            cfeature.LAKES,
            linewidth=1,
            linestyle="-",
            edgecolor="dimgrey",
            facecolor="none",
        )
        # Actually make the plot
        if style == None:
            c = ax.pcolormesh(
                data[lon_dim],
                data[lat_dim],
                np.mean(data["precipitation"][0, :, valid_time_idx, :, :], axis=0)
                * plot_norm,
                norm=colors.LogNorm(*value_range_precip),
                cmap="YlGnBu",
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(
                data[lon_dim],
                data[lat_dim],
                np.mean(data["precipitation"][0, :, valid_time_idx, :, :], axis=0)
                * plot_norm,
                colors=plot_colours,
                levels=plot_levels * plot_norm,
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(
                c, fraction=0.04, ticks=plot_levels * plot_norm
            )  # Add a colorbar with a nice size
        if region != COUNTRY_NAMES[0]:
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar
        ax.set_title(f"Ensemble mean", size=14)  # This plot's title

        ax = axs[1]  # Second plot (right)
        ax.add_feature(
            cfeature.COASTLINE, linewidth=1
        )  # Draw some features to see where we are
        ax.add_feature(shape_feature)  # The borders
        ax.add_feature(
            cfeature.LAKES,
            linewidth=1,
            linestyle="-",
            edgecolor="dimgrey",
            facecolor="none",
        )
        # Actually make the plot
        if style == None:
            c = ax.pcolormesh(
                data[lon_dim],
                data[lat_dim],
                np.std(
                    data["precipitation"][0, :, valid_time_idx, :, :], axis=0, ddof=1
                )
                * plot_norm,
                norm=colors.LogNorm(*value_range_precip),
                cmap="YlGnBu",
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(
                data[lon_dim],
                data[lat_dim],
                np.std(
                    data["precipitation"][0, :, valid_time_idx, :, :], axis=0, ddof=1
                )
                * plot_norm,
                colors=plot_colours,
                levels=plot_levels * plot_norm,
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(
                c, fraction=0.04, ticks=plot_levels * plot_norm
            )  # Add a colorbar with a nice size
        if region != COUNTRY_NAMES[0]:
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar
        ax.set_title(f"Ensemble standard deviation", size=14)  # This plot's title

        fig.suptitle(
            f"Jurre Brishti cGAN forecast: Valid {valid_time} - {valid_time + timedelta(hours=6)}"
        )  # Overall title
        plt.tight_layout()  # Looks nicer
        plt.show()  # Finally draw the plot


# Plot all ensemble members in the cGAN forecast data at a specified valid time.
# Arguments
#   data                  - An xarray DataSet containing the cGAN rainfall forecasts.
#   valid_time_start_hour - The hour the valid time starts at. Can either be 6, 12, 18 or 0.
#   style=None            - Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
#   plot_units='mm/h'     - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   region='ICPAC'        - can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                           'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
def plot_GAN_ensemble(
    data: xr.Dataset,
    valid_time_start_hour: int,
    lon_dim: str | None = "longitude",
    lat_dim: str | None = "latitude",
    style: str | None = COLOR_SCHEMES[0],
    plot_units: str | None = "mm/h",
    region: str | None = COUNTRY_NAMES[0],
):

    # Get the units to use for plotting
    plot_norm = get_plot_normalisation(plot_units)

    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15 * plot_norm)

    # Use a style other than the default
    if style is not None:
        plot_levels, plot_colours = get_contour_levels(style)

    # Load the border shapefile
    reader = get_shape_boundary(shape_name=region)
    shape_feature = ShapelyFeature(
        reader.geometries(), ccrs.PlateCarree(), facecolor="none"
    )

    # Get the extent of the region that we are looking at
    if region != COUNTRY_NAMES[0]:
        region_extent = get_region_extent(region, border_size=0.5)

    # Change valid_time_start_hour into the valid_time_idx
    if valid_time_start_hour == 6:
        valid_time_idx = 0
    elif valid_time_start_hour == 12:
        valid_time_idx = 1
    elif valid_time_start_hour == 18:
        valid_time_idx = 2
    elif valid_time_start_hour == 0:
        valid_time_idx = 3
    else:
        print("ERROR: valid_time_start_hour must be 6, 12, 18 or 0.")
        return

    # Convert the forecast valid time to a datetime.datetime format
    valid_time = datetime64_to_datetime(
        data["fcst_valid_time"][0, valid_time_idx].values
    )

    # Define the figure and each axis for the rows and columns
    fig, axs = plt.subplots(
        nrows=10,
        ncols=5,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(10, 25),
        layout="constrained",
    )

    # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs = axs.flatten()

    # For each ensemble member
    for ax_idx in range(data["member"].size):

        ax = axs[ax_idx]  # First plot (left)
        ax.add_feature(
            cfeature.COASTLINE, linewidth=1
        )  # Draw some features to see where we are
        ax.add_feature(shape_feature)  # The borders
        ax.add_feature(
            cfeature.LAKES,
            linewidth=1,
            linestyle="-",
            edgecolor="dimgrey",
            facecolor="none",
        )
        # Actually make the plot
        if style == None:
            c = ax.pcolormesh(
                data[lon_dim],
                data[lat_dim],
                data["precipitation"][0, ax_idx, valid_time_idx, :, :] * plot_norm,
                norm=colors.LogNorm(*value_range_precip),
                cmap="YlGnBu",
                transform=ccrs.PlateCarree(),
            )
        else:
            c = ax.contourf(
                data[lon_dim],
                data[lat_dim],
                data["precipitation"][0, ax_idx, valid_time_idx, :, :] * plot_norm,
                colors=plot_colours,
                levels=plot_levels * plot_norm,
                transform=ccrs.PlateCarree(),
            )
        if region != COUNTRY_NAMES[0]:
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        ax.set_title(f"{ax_idx+1}", size=14)  # This plot's title

    # Add a final colorbar with a nice size
    if style == None:
        cb = fig.colorbar(c, ax=axs, location="bottom", shrink=0.4, pad=0.01)
    else:
        cb = fig.colorbar(
            c,
            ax=axs,
            location="bottom",
            shrink=0.4,
            pad=0.01,
            ticks=plot_levels * plot_norm,
        )
    cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar

    fig.suptitle(
        f"Jurre Brishti cGAN ensemble: Valid {valid_time} - {valid_time + timedelta(hours=6)}"
    )  # Overall title
    plt.show()  # Finally draw the plot
