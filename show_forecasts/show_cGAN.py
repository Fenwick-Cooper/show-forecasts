# Load and plot cGAN forecast data.

# To do:
#   Change "mm h**-1" to "mm/h" in the data.
#   Store and extract the model name from the data.
#   Add the initialisation time to the title.
from os import getenv
from typing import Literal
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
import matplotlib.pyplot as plt
from matplotlib import colors  # For consistency with Harris et. al 2022
from datetime import datetime, date
from pathlib import Path
import xarray as xr
from show_forecasts.data_utils import (
    get_region_extent,
    get_shape_boundary,
    get_plot_normalisation,
    get_contour_levels,
    datetime64_to_datetime,
    get_threshold_plot_colours,
)
from show_forecasts.constants import (
    COUNTRY_NAMES,
    COLOR_SCHEMES,
    GAN_THRESHOLD_PLOT_LEVEL_NAMES,
    GAN_THRESHOLD_PLOT_LEVELS,
)


# Load a 24 hour mean forecast at a lead time of 30 to 54 hours
# To be clear about the dates use print_forecast_info(forecast_init_date)
# Arguments
#    data_dir                 - Directory where the data is stored.
#    mask_region (optional)   - region to be used for plotting. defaults to East Africa
#                              valid options are 'East Africa', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                              'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#    init_date                - A datetime.datetime corresponding to when the forecast was initialised.
#    init_time                - A two digits string to denote forecast initialization time. Valid values are 00, 06, 12 and 18
# cgan_ui_fs (optional)       - instruction on whether to use new cgan ui filesystem structure. defaults to false
# Returns
#    An xarray DataSet containing the cGAN rainfall forecasts.
def load_GAN_forecast(
    model: str,
    data_dir: str,
    init_date: datetime | date,
    init_time: Literal["00", "06", "12", "18"] | None = "00",
    mask_region: str | None = COUNTRY_NAMES[0],
    cgan_ui_fs: bool | None = False,
) -> xr.Dataset:
    init_time = "00" if init_time not in ["06", "12", "18"] else init_time
    if cgan_ui_fs:
        mask_region = mask_region if mask_region is not None else COUNTRY_NAMES[0]
        fcst_filename = (
            f"{mask_region.lower().replace(' ','_')}-{model.replace('-','_')}-"
            + f"{init_date.year}{init_date.month:02}{init_date.day:02}_{init_time}Z.nc"
        )
        file_path = (
            Path(data_dir)
            / model
            / mask_region
            / str(init_date.year)
            / str(init_date.month).rjust(2, "0")
            / fcst_filename
        )
    else:
        file_path = f"{data_dir}/GAN_{init_date.year}{init_date.month:02}{init_date.day:02}_{init_time}Z.nc"
    data = xr.open_dataset(file_path)
    return data


# Sorts the data along the ensemble axis. Used in percentile plots.
# Arguments
#   data                    - An xarray DataSet containing the cGAN rainfall forecasts.
# Returns
#   data_sorted             - A xarray DataSet containing the cGAN rainfall forecasts sorted along the
#                             ensemble axis.
def sort_along_ensemble_axis(data):

    # Create a copy of the Dataset to hold the sorted data
    data_sorted = data.copy()

    # Find the axis to sort along
    member_axis = data["precipitation"].get_axis_num("member")

    # Sort values along the chosen axis
    precip_sorted = np.sort(data["precipitation"], axis=member_axis)

    # Copy the sorted data into the data_sorted Dataset
    data_sorted["precipitation"].values = precip_sorted

    return data_sorted


# Plot the ensemble mean and ensemble standard deviation of the cGAN forecast data
# at each valid time.
# Arguments
#   data                         - An xarray DataSet containing the cGAN rainfall forecasts.
#   model                        - name of cGAN model. One of jurre-brishti-ens of mvua-kubwa-ens
#   accumulation_time='06h'      - Can be '06h', or '24h'
#   valid_time_start_hour='all' - The hour the valid time starts at.
#   style=None                  - Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
#   plot_units='mm/h'           - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   region='East Africa'              - can be 'East Africa', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi',
#                                 'Djibouti', 'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#   file_name=None              - If a file name, ending in '.png', '.jpg' or '.pdf' is
#                                 specified, the plot is saved in that format. If
#                                 valid_time_start_hour = 'all', the hour is appended to the
#                                 file name.
def plot_GAN_forecast(
    data: xr.Dataset,
    model: Literal["jurre-brishti-ens", "mvua-kubwa-ens"] | None = "jurre-brishti-ens",
    accumulation_time: str | None = "06h",
    valid_time_start_hour: str | None = "all",
    lon_dim: str | None = "longitude",
    lat_dim: str | None = "latitude",
    style: str | None = COLOR_SCHEMES[0],
    plot_units: str | None = "mm/h",
    region: str | None = COUNTRY_NAMES[0],
    file_name: str | None = None,
    show_plot: bool | None = True,
):

    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)

    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15 * plot_norm)

    # Use a style other than the default
    if style is not None:
        plot_levels, plot_colours = get_contour_levels(style)

    # Load EA region border shapefile
    reader = get_shape_boundary()
    shape_feature = ShapelyFeature(
        reader.geometries(), ccrs.PlateCarree(), facecolor="none"
    )

    # Get the extent of the region that we are looking at
    if region != COUNTRY_NAMES[0] and region is not None:
        # load region of interest boundary layer
        reader = get_shape_boundary(shape_name=region)
        region_feature = ShapelyFeature(
            reader.geometries(), ccrs.PlateCarree(), facecolor="none"
        )
        # Get extents of the region we are looking at
        region_extent = get_region_extent(region, border_size=0.5)

    forecast_valid_times = {
        "jurre-brishti-ens": ["30h", "36h", "42h", "48h"],
        "mvua-kubwa-ens": ["06h", "30h", "54h", "78h", "102h", "126h", "150h"],
    }
    # set default forecast valid start time
    valid_time_idx_list = [0]

    if accumulation_time == "06h":

        # Change valid_time_start_hour into the valid_time_idx
        if valid_time_start_hour == "all":
            valid_time_idx_list = range(len(data["valid_time"]))
        elif model in forecast_valid_times.keys():
            if valid_time_start_hour in forecast_valid_times[model]:
                valid_time_idx_list = [
                    forecast_valid_times[model].index(valid_time_start_hour)
                ]
            else:
                print(
                    f"ERROR: valid_time_start_hour must be one of {', '.join(forecast_valid_times[model])} or 'all' for {model} model."
                )
        else:
            print(
                f"ERROR: model name must be one of jurre-brishti-ens or mvua-kubwa-ens. found invalid {model} model name"
            )

    elif accumulation_time == "24h":

        if (valid_time_start_hour != "30h") and (valid_time_start_hour != "all"):
            print(
                "ERROR: valid_time_start_hour must be 30h when accumulation_time is '24h'."
            )
            valid_time_start_hour = "30h"

    else:
        print("ERROR: accumulation_time must be either '06h' or '24h'.")

    # There are plots for each valid time
    for valid_time_idx in valid_time_idx_list:

        # Convert the forecast initialization time to a datetime.datetime format
        fcst_init_time = datetime64_to_datetime(data["time"][0].values)

        # Convert the forecast valid time to a datetime.datetime format
        fcst_valid_time = datetime64_to_datetime(
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
        ax.add_feature(
            cfeature.LAKES,
            linewidth=1,
            linestyle="-",
            edgecolor="dimgrey",
            facecolor="none",
        )
        ax.add_feature(shape_feature)  # EA region borders
        if region != COUNTRY_NAMES[0] and region is not None:
            ax.add_feature(region_feature, linestyle=":")
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        # Either plot 6h data or 24h data
        if accumulation_time == "06h":
            data_to_plot = (
                np.mean(data["precipitation"][0, :, valid_time_idx, :, :], axis=0)
                * plot_norm
            )
        elif accumulation_time == "24h":
            data_to_plot = (
                np.mean(data["precipitation"][0, :, :, :, :], axis=(0, 1)) * plot_norm
            )
        # Actually make the plot
        if style is None:
            c = ax.pcolormesh(
                data[lon_dim],
                data[lat_dim],
                data_to_plot,
                norm=colors.LogNorm(*value_range_precip),
                cmap="YlGnBu",
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(
                data[lon_dim],
                data[lat_dim],
                data_to_plot,
                colors=plot_colours,
                levels=plot_levels * plot_norm,
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
            cb_labels = np.round(plot_levels * plot_norm, 1).astype(str).tolist()
            cb_labels[-1] = ""  # Remove the final label
            cb.set_ticks(ticks=plot_levels * plot_norm, labels=cb_labels)
        cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar
        ax.set_title("Ensemble mean", size=14)  # This plot's title

        ax = axs[1]  # Second plot (right)
        ax.add_feature(
            cfeature.COASTLINE, linewidth=1
        )  # Draw some features to see where we are
        ax.add_feature(
            cfeature.LAKES,
            linewidth=1,
            linestyle="-",
            edgecolor="dimgrey",
            facecolor="none",
        )
        ax.add_feature(shape_feature)  # EA region borders
        if region != COUNTRY_NAMES[0] and region is not None:
            ax.add_feature(region_feature, linestyle=":")
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        # Either plot 6h data or 24h data
        if accumulation_time == "06h":
            data_to_plot = (
                np.std(
                    data["precipitation"][0, :, valid_time_idx, :, :], axis=0, ddof=1
                )
                * plot_norm
            )
        elif accumulation_time == "24h":
            data_to_plot = (
                np.sqrt(
                    np.mean(
                        np.var(data["precipitation"][0, :, :, :, :], axis=0, ddof=1),
                        axis=0,
                    )
                )
                * plot_norm
            )
        # Actually make the plot
        if style is None:
            c = ax.pcolormesh(
                data[lon_dim],
                data[lat_dim],
                data_to_plot,
                norm=colors.LogNorm(*value_range_precip),
                cmap="YlGnBu",
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(
                data[lon_dim],
                data[lat_dim],
                data_to_plot,
                colors=plot_colours,
                levels=plot_levels * plot_norm,
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
            cb_labels = np.round(plot_levels * plot_norm, 1).astype(str).tolist()
            cb_labels[-1] = ""  # Remove the final label
            cb.set_ticks(ticks=plot_levels * plot_norm, labels=cb_labels)
        cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar
        ax.set_title("Ensemble standard deviation", size=14)  # This plot's title

        fig.suptitle(
            f"{model.replace('-', '  ').replace('ens','').title()} cGAN forecast: Valid {fcst_init_time.strftime('%Y-%m-%d %H:00')} to {fcst_valid_time.strftime('%Y-%m-%d %H:00')} {getenv('DEFAULT_TIMEZONE', 'UTC')}"
        )  # Overall title
        plt.tight_layout()  # Looks nicer

        # Save the plot
        if file_name is not None:
            if file_name[-4:] in [".png", ".jpg", ".pdf"]:
                # If we are making more than one plot
                if valid_time_start_hour == "all":
                    # Append the hour to the file name
                    save_file_name = f"{file_name[:-4]}_{fcst_init_time.hour:02d}_{fcst_valid_time.hour:02d}{file_name[-4:]}"
                else:  # We are making only one plot
                    save_file_name = file_name  # Use the exact file name specified
                plt.savefig(save_file_name, format=file_name[-3:], bbox_inches="tight")
            else:
                print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")

        if show_plot:
            plt.show()  # Finally draw the plot


# Plot all ensemble members in the cGAN forecast data at a specified valid time.
# Arguments
#   data                  - An xarray DataSet containing the cGAN rainfall forecasts.
#   model                        - name of cGAN model. One of jurre-brishti-ens of mvua-kubwa-ens
#   valid_time_start_hour - The hour the valid time starts at.
#   style=None            - Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
#   plot_units='mm/h'     - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   region='East Africa'  - can be 'East Africa', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                           'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#   max_num_plots=50      - Maximum number of ensemble members to plot.
#   file_name=None        - If a file name, ending in '.png', '.jpg' or '.pdf' is specified, the
#                            plot is saved in that format.
def plot_GAN_ensemble(
    data: xr.Dataset,
    model: Literal["jurre-brishti-ens", "mvua-kubwa-ens"] | None = "jurre-brishti-ens",
    valid_time_start_hour: str | None = "30",
    lon_dim: str | None = "longitude",
    lat_dim: str | None = "latitude",
    style: str | None = COLOR_SCHEMES[0],
    plot_units: str | None = "mm/h",
    region: str | None = COUNTRY_NAMES[0],
    max_num_plots: int | None = 50,
    file_name: str | None = None,
    show_plot: bool | None = True,
):
    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)

    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15 * plot_norm)

    # Use a style other than the default
    if style is not None:
        plot_levels, plot_colours = get_contour_levels(style)

    # Load EA region border shapefile
    reader = get_shape_boundary()
    shape_feature = ShapelyFeature(
        reader.geometries(), ccrs.PlateCarree(), facecolor="none"
    )
    # Get the extent of the region that we are looking at
    if region != COUNTRY_NAMES[0] and region is not None:
        # load selected area of interest boundary
        reader = get_shape_boundary(shape_name=region)
        region_feature = ShapelyFeature(
            reader.geometries(), ccrs.PlateCarree(), facecolor="none"
        )
        # Get extents of the region we are looking at
        region_extent = get_region_extent(region, border_size=0.5)

    # Change valid_time_start_hour into the valid_time_idx
    forecast_valid_times = {
        "jurre-brishti-ens": ["30h", "36h", "42h", "48h"],
        "mvua-kubwa-ens": ["06h", "30h", "54h", "78h", "102h", "126h", "150h"],
    }
    # set default forecast valid start time
    valid_time_idx = 0

    # Change valid_time_start_hour into the valid_time_idx
    if model in forecast_valid_times.keys():
        if valid_time_start_hour in forecast_valid_times[model]:
            valid_time_idx = forecast_valid_times[model].index(valid_time_start_hour)
        else:
            print(
                f"ERROR: valid_time_start_hour must be one of {', '.join(forecast_valid_times[model])} or 'all' for {model} model."
            )
    else:
        print(
            f"ERROR: model name must be one of jurre-brishti-ens or mvua-kubwa-ens. found invalid {model} model name"
        )

    # Convert the forecast initialization time to a datetime.datetime format
    fcst_init_time = datetime64_to_datetime(data["time"][0].values)

    # Convert the forecast valid time to a datetime.datetime format
    fcst_valid_time = datetime64_to_datetime(
        data["fcst_valid_time"][0, valid_time_idx].values
    )

    # How many plots will we make
    num_plots = np.min([max_num_plots, data["member"].size])
    num_rows = int(np.ceil(num_plots / 5))

    # Define the figure and each axis for the rows and columns
    fig, axs = plt.subplots(
        nrows=num_rows,
        ncols=5,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(10, 2.8 * num_rows + 0.9),
        layout="constrained",
    )

    # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs = axs.flatten()

    # Don't show axes without plots
    for ax_idx in range(num_plots, num_rows * 5):
        axs[ax_idx].set_axis_off()

    # Don't show axes without plots
    for ax_idx in range(num_plots, num_rows * 5):
        axs[ax_idx].set_axis_off()

    # For each ensemble member
    for ax_idx in range(num_plots):

        ax = axs[ax_idx]  # First plot (left)
        ax.add_feature(
            cfeature.COASTLINE, linewidth=1
        )  # Draw some features to see where we are
        ax.add_feature(
            cfeature.LAKES,
            linewidth=1,
            linestyle="-",
            edgecolor="dimgrey",
            facecolor="none",
        )
        ax.add_feature(shape_feature)  # The EA region borders
        if region != COUNTRY_NAMES[0] and region is not None:
            ax.add_feature(region_feature, linestyle=":")
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        # Actually make the plot
        if style is None:
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
        ax.set_title(f"{ax_idx+1}", size=14)  # This plot's title

    # Add a final colorbar with a nice size
    cb = fig.colorbar(c, ax=axs, location="bottom", shrink=0.6, pad=0.01)
    if style is not None:
        cb_labels = np.round(plot_levels * plot_norm, 1).astype(str).tolist()
        cb_labels[-1] = ""  # Remove the final label
        cb.set_ticks(ticks=plot_levels * plot_norm, labels=cb_labels)
    cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar

    fig.suptitle(
        f"{model.replace('-', '  ').replace('ens','').title()} cGAN ensemble: Valid {fcst_init_time.strftime('%Y-%m-%d %H:00')} to {fcst_valid_time.strftime('%Y-%m-%d %H:00')} {getenv('DEFAULT_TIMEZONE', 'UTC')}"
    )  # Overall title

    # Save the plot
    if file_name is not None:
        if file_name[-4:] in [".png", ".jpg", ".pdf"]:
            plt.savefig(
                f"{file_name[:-4]}_{fcst_init_time.hour:02d}{file_name[-4:]}",
                format=file_name[-3:],
                bbox_inches="tight",
            )
        else:
            print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")

    if show_plot:
        plt.show()  # Finally draw the plot


# Plot the chance of rainfall at a rate above a specified threshold.
# Arguments
#   data                    - An xarray DataSet containing the cGAN rainfall forecasts.
#   model                        - name of cGAN model. One of jurre-brishti-ens of mvua-kubwa-ens
#   threshold=2             - We'll plot the chance of rainfall above this threshold rate. The
#                             default is 2 mm/h. The units of threshold is set by plot_units.
#   plot_units='mm/h'       - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   valid_time_start_hour=30 - The hour the valid time starts at. Can either be 6, 12, 18 or 0 UTC,
#                             or specify 'all' to make all four plots.
#   style=None              - Options: 'ICPAC', 'KMD', 'EMI'
#   show_percentages=False  - Either shows a description (False) or the percentage (True) of
#                             the chance of exceeding the threshold.
#   region='East Africa'    - Can be 'East Africa', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                             'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#   file_name=None          - If a file name, ending in '.png', '.jpg' or '.pdf' is specified, the
#                             plot is saved in that format.
def plot_GAN_threshold_chance(
    data: xr.Dataset,
    model: Literal["jurre-brishti-ens", "mvua-kubwa-ens"] | None = "jurre-brishti-ens",
    valid_time_start_hour: str | None = "all",
    style: str | None = "ICPAC",
    threshold: str | None = 2,
    plot_units: str | None = "mm/h",
    show_percentages: bool | None = False,
    region: str | None = COUNTRY_NAMES[0],
    file_name: str | None = None,
    show_plot: bool | None = True,
):

    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)

    # Normalise the threshold to the chosen units
    threshold /= plot_norm
    plot_level_percentages = []
    for i in range(len(GAN_THRESHOLD_PLOT_LEVELS)):
        plot_level_percentages.append(f"{GAN_THRESHOLD_PLOT_LEVELS[i]}%")

    # Load the border shapefile
    reader = get_shape_boundary()
    borders_feature = ShapelyFeature(
        reader.geometries(), ccrs.PlateCarree(), facecolor="none"
    )

    if region != COUNTRY_NAMES[0] and region is not None:

        # Load the regions shapefile
        reader = get_shape_boundary(shape_name=region)
        regions_feature = ShapelyFeature(
            reader.geometries(), ccrs.PlateCarree(), facecolor="none"
        )

        # Get the extent of the region that we are looking at
        region_extent = get_region_extent(region, border_size=0.5)

    forecast_valid_times = {
        "jurre-brishti-ens": ["30h", "36h", "42h", "48h"],
        "mvua-kubwa-ens": ["06h", "30h", "54h", "78h", "102h", "126h", "150h"],
    }
    # set default forecast valid start time
    valid_time_idx_list = [0]

    # Change valid_time_start_hour into the valid_time_idx
    if valid_time_start_hour == "all":
        valid_time_idx_list = range(len(data["valid_time"]))
    elif model in forecast_valid_times.keys():
        if valid_time_start_hour in forecast_valid_times[model]:
            valid_time_idx_list = [
                forecast_valid_times[model].index(valid_time_start_hour)
            ]
        else:
            print(
                f"ERROR: valid_time_start_hour must be one of {', '.join(forecast_valid_times[model])} or 'all' for {model} model."
            )
    else:
        print(
            f"ERROR: model name must be one of jurre-brishti-ens or mvua-kubwa-ens. found invalid {model} model name"
        )

    if len(valid_time_idx_list) == 1:

        # Define the figure and axes
        fig, axs = plt.subplots(
            nrows=1,
            ncols=1,
            subplot_kw={"projection": ccrs.PlateCarree()},
            figsize=(5, 5),
        )

        # axs is a `GeoAxes`. Make it into a 1-D array
        axs = [axs]

    else:
        figsize = (8, 8) if len(valid_time_idx_list) == 4 else (12, 18)
        # Define the figure and axes
        fig, axs = plt.subplots(
            nrows=int(np.ceil(len(valid_time_idx_list) / 2)),
            ncols=2,
            subplot_kw={"projection": ccrs.PlateCarree()},
            figsize=figsize,
        )

        # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
        axs = axs.flatten()

    # Define the plot colours
    plot_colours = get_threshold_plot_colours(style)

    # Define the plot colours
    plot_colours = get_threshold_plot_colours(style)

    # Convert the forecast valid time to a datetime.datetime format
    fcst_init_time = datetime64_to_datetime(data["time"][0].values)

    # There are plots for each valid time
    for idx, valid_time_idx in enumerate(valid_time_idx_list):

        # Convert the forecast valid time to a datetime.datetime format
        fcst_valid_time = datetime64_to_datetime(
            data["fcst_valid_time"][0, valid_time_idx].values
        )

        # Keep the first valid time for the plot title
        # TODO: remove this line
        # if idx == 0:
        #     first_valid_time = fcst_valid_time

        # The percentage of ensemble members that exceed the threshold
        plot_data = (
            np.sum(
                data["precipitation"][0, :, valid_time_idx, :, :] > threshold, axis=0
            )
            * 100
            / len(data["member"])
        )

        ax = axs[idx]
        ax.gridlines()
        ax.set_facecolor("white")  # For consistency with Harris et. al 2022
        ax.add_feature(cfeature.COASTLINE, linewidth=1)
        if region != COUNTRY_NAMES[0] and region is not None:
            ax.add_feature(regions_feature, linestyle=":")
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        ax.add_feature(borders_feature)  # The borders
        ax.add_feature(
            cfeature.LAKES,
            linewidth=1,
            linestyle="-",
            edgecolor="dimgrey",
            facecolor="none",
        )
        c = ax.contourf(
            data["longitude"],
            data["latitude"],
            plot_data,
            levels=GAN_THRESHOLD_PLOT_LEVELS,
            transform=ccrs.PlateCarree(),
            colors=plot_colours,
        )
        ax.set_title(
            f"{fcst_init_time.strftime('%Y-%m-%d %H:00')} - {fcst_valid_time.strftime('%Y-%m-%d %H:00')} {getenv('DEFAULT_TIMEZONE', 'UTC')}",
            size=14,
        )
        cb = plt.colorbar(c, fraction=0.04)
        # cb.ax.tick_params(labelsize=18)
        if show_percentages:
            # cb.set_label(f'% chance',size=18)
            cb.set_ticks(ticks=GAN_THRESHOLD_PLOT_LEVELS, labels=plot_level_percentages)
        else:
            cb.set_ticks(
                ticks=GAN_THRESHOLD_PLOT_LEVELS, labels=GAN_THRESHOLD_PLOT_LEVEL_NAMES
            )

    title_string = f"""{model.replace('-', '  ').replace('ens','').title()} Threshold Chance: Valid {fcst_init_time.strftime('%Y-%m-%d %H:00')} to {fcst_valid_time.strftime('%Y-%m-%d %H:00')}
    Chance of rainfall above {threshold*plot_norm:.1f} {plot_units}."""

    fig.suptitle(title_string)  # Overall title
    plt.tight_layout()  # Looks nicer

    # Save the plot
    if file_name is not None:
        if file_name[-4:] in [".png", ".jpg", ".pdf"]:
            plt.savefig(file_name, format=file_name[-3:], bbox_inches="tight")
        else:
            print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")

    if show_plot:
        plt.show()
