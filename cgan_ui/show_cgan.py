# Load and plot cGAN forecast data.

# To do:
#   Change "mm h**-1" to "mm/h" in the data.
#   Store the model name in the data.
#   Combine with the show_forecasts.py script used for plotting IFS open data.

import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import colors  # For consistency with Harris et. al 2022
from datetime import datetime, timedelta
import xarray as xr


# Convert numpy.datetime64 to datetime
# Arguments
#    datetime64 - A single date stored in a datetime64 object.
# Returns
#    A date as a datetime.datetime object.
def datetime64_to_datetime(datetime64):
    unix_epoch = np.datetime64(0, "s")
    one_second = np.timedelta64(1, "s")
    seconds_since_epoch = (datetime64 - unix_epoch) / one_second
    return datetime.utcfromtimestamp(seconds_since_epoch)


# Returns the contour levels and colours used for different styles of plot
# Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
# XXX Incorrect levels for ICPAC Heavy
# XXX Incorrect levels for EMI Heavy
def get_contour_levels(style):

    if style == "ICPAC_heavy":

        # ICPAC Heavy
        # 00008b rgb(0, 0, 139)      Extremely Heavy Rainfall
        # 1874cd rgb(24, 116, 205)   Very Heavy Rainfall
        # afeeee rgb(175, 238, 238)  Heavy Rainfall
        # ffffff rgb(255, 255, 255)  Less than heavy

        # XXX Incorrect levels for ICPAC Heavy
        plot_levels = np.array([0, 50, 100, 200, 400]) / (
            7 * 24
        )  # Convert from mm/week to mm/h
        plot_colours = [
            [1.0, 1.0, 1.0, 1.0],
            [0.68627451, 0.93333333, 0.93333333, 1.0],
            [0.09411765, 0.45490196, 0.80392157, 1.0],
            [0.0, 0.0, 0.54509804, 1.0],
        ]

    elif style == "ICPAC":

        # ICPAC
        # 228b22 rgb(34, 139, 34)    Above 200 mm/day
        # 6cd403 rgb(108, 212, 3)    100-200 mm/day
        # 00fe00 rgb(0, 254, 0)      50-100 mm/day
        # caff70 rgb(202, 255, 112)  30-50 mm/day
        # ffff00 rgb(255, 255, 0)    10-30 mm/day
        # ffa500 rgb(255, 165, 0)    1-10 mm/day
        # d9d9d9 rgb(217, 217, 217)  Less than 1 mm/day

        plot_levels = np.array([0, 1, 10, 30, 50, 100, 200, 800]) / (
            7 * 24
        )  # Convert from mm/week to mm/h
        plot_colours = [
            [0.85098039, 0.85098039, 0.85098039, 1.0],
            [1.0, 0.64705882, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
            [0.79215686, 1.0, 0.43921569, 1.0],
            [0.0, 0.99607843, 0.0, 1.0],
            [0.42352941, 0.83137255, 0.01176471, 1.0],
            [0.13333333, 0.54509804, 0.13333333, 1.0],
        ]

    elif style == "KMD":

        # KMD
        # 910000 rgb(145, 0, 0)      Above 100 mm/day
        # ee0005 rgb(238, 0, 5)      80-100 mm/day
        # ee6f01 rgb(238, 111, 1)    70-80 mm/day
        # faa700 rgb(250, 167, 0)    60-70 mm/day
        # ffc800 rgb(255, 200, 0)    50-60 mm/day
        # 2481c5 rgb(36, 129, 197)   40-50 mm/day
        # 45a7e6 rgb(69, 167, 230)   35-40 mm/day
        # 8fc2e5 rgb(143, 194, 229)  30-35 mm/day
        # b6d9fc rgb(182, 217, 252)  25-30 mm/day
        # ccebfd rgb(204, 235, 253)  20-25 mm/day
        # 1a8f1a rgb(26, 143, 26)    15-20 mm/day
        # 2cba28 rgb(44, 186, 40)    10-15 mm/day
        # 39d904 rgb(57, 217, 4)     5-10 mm/day
        # 3ef600 rgb(62, 246, 0)     2-5 mm/day
        # ffffff rgb(255, 255, 255)  Less than 1 mm/day

        plot_levels = (
            np.array([0, 2, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 200])
            / 24
        )  # Convert from mm/day to mm/h
        plot_colours = [
            [1.0, 1.0, 1.0, 1.0],
            [0.24313725, 0.96470588, 0.0, 1.0],
            [0.22352941, 0.85098039, 0.01568627, 1.0],
            [0.17254902, 0.72941176, 0.15686275, 1.0],
            [0.10196078, 0.56078431, 0.10196078, 1.0],
            [0.8, 0.92156863, 0.99215686, 1.0],
            [0.71372549, 0.85098039, 0.98823529, 1.0],
            [0.56078431, 0.76078431, 0.89803922, 1.0],
            [0.27058824, 0.65490196, 0.90196078, 1.0],
            [0.14117647, 0.50588235, 0.77254902, 1.0],
            [1.0, 0.78431373, 0.0, 1.0],
            [0.98039216, 0.65490196, 0.0, 1.0],
            [0.93333333, 0.43529412, 0.00392157, 1.0],
            [0.93333333, 0.0, 0.01960784, 1.0],
            [0.56862745, 0.0, 0.0, 1.0],
        ]

    elif style == "EMI_heavy":

        # EMI Heavy
        # 00008b rgb(0, 0, 139)      Extremely Heavy Rainfall
        # 1874cd rgb(24, 116, 205)   Very Heavy Rainfall
        # afeeee rgb(175, 238, 238)  Heavy Rainfall
        # ffffff rgb(255, 255, 255)  Less than heavy

        # XXX Incorrect levels for EMI Heavy
        plot_levels = np.array([0, 50, 100, 200, 400]) / (
            7 * 24
        )  # Convert from mm/week to mm/h
        plot_colours = [
            [1.0, 1.0, 1.0, 1.0],
            [0.68627451, 0.93333333, 0.93333333, 1.0],
            [0.09411765, 0.45490196, 0.80392157, 1.0],
            [0.0, 0.0, 0.54509804, 1.0],
        ]

    elif style == "EMI":

        # EMI
        # 228b22 rgb(34, 139, 34)    Above 100 mm/day
        # 00ff00 rgb(0, 255, 0)      75-100 mm/day
        # 66cd00 rgb(102, 205, 0)    50-75 mm/day
        # 7fff00 rgb(127, 255, 0)    30-50 mm/day
        # a2cd5a rgb(162, 205, 90)   20-30 mm/day
        # caff70 rgb(202, 255, 112)  10-20 mm/day
        # ffff00 rgb(255, 255, 0)    5-10 mm/day
        # ffa500 rgb(255, 165, 0)    1-5 mm/day
        # d9d9d9 rgb(217, 217, 217)  Less than 1 mm/day

        plot_levels = (
            np.array([0, 1, 5, 10, 20, 30, 50, 75, 100, 200]) / 24
        )  # Convert from mm/day to mm/h
        plot_colours = [
            [0.85098039, 0.85098039, 0.85098039, 1.0],
            [1.0, 0.64705882, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
            [0.79215686, 1.0, 0.43921569, 1.0],
            [0.63529412, 0.80392157, 0.35294118, 1.0],
            [0.49803922, 1.0, 0.0, 1.0],
            [0.4, 0.80392157, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [0.13333333, 0.54509804, 0.13333333, 1.0],
        ]

    else:
        print(f"ERROR: Unknown style {style}")

    return plot_levels, plot_colours


# Load a 24 hour mean forecast at a lead time of 30 to 54 hours
# To be clear about the dates use print_forecast_info(forecast_init_date)
# Arguments
#    forecast_init_date   - A datetime.datetime corresponding to when the forecast was initialised.
#    data_dir             - Directory where the data is stored.
# Returns
#    An xarray DataSet containing the cGAN rainfall forecasts.
def load_GAN_forecast(forecast_init_date, data_dir):
    d = forecast_init_date  # Shorthand
    file_name = f"{data_dir}/GAN_{d.year}{d.month:02}{d.day:02}.nc"
    data = xr.open_dataset(file_name)
    return data


# Plot the ensemble mean and ensemble standard deviation of the cGAN forecast data
# at each valid time.
# Arguments
#   data              - An xarray DataSet containing the cGAN rainfall forecasts.
#   style=None        - Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
#   plot_units='mm/h' - Can be 'mm/h' (default), 'mm/day' or 'mm/week'
def plot_GAN_forecast(data, style=None, plot_units="mm/h"):

    if plot_units == "mm/h":
        plot_norm = 1
    elif plot_units == "mm/day":
        plot_norm = 24
    elif plot_units == "mm/week":
        plot_norm = 7 * 24
    else:
        print(f"ERROR: Unknown plot units {plot_units}")
        print(f"       Options are 'mm/h', 'mm/day', 'mm/week'.")
        return

    # Use a style other than the default
    if style is not None:
        plot_levels, plot_colours = get_contour_levels(style)

    # There are plots for each valid time
    for valid_time_idx in range(len(data["valid_time"])):

        # Convert the forecast valid time to a datetime.datetime format
        valid_time = datetime64_to_datetime(
            data["fcst_valid_time"][0, valid_time_idx].values
        )

        # To be consistent with the Harris et. al paper.
        value_range_precip = (0.1, 15 * plot_norm)

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
        ax.add_feature(cfeature.BORDERS, linewidth=1)
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
                data["longitude"],
                data["latitude"],
                np.mean(data["precipitation"][0, :, valid_time_idx, :, :], axis=0)
                * plot_norm,
                norm=colors.LogNorm(*value_range_precip),
                cmap="YlGnBu",
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(
                data["longitude"],
                data["latitude"],
                np.mean(data["precipitation"][0, :, valid_time_idx, :, :], axis=0)
                * plot_norm,
                colors=plot_colours,
                levels=plot_levels * plot_norm,
                transform=ccrs.PlateCarree(),
            )
            cb = plt.colorbar(
                c, fraction=0.04, ticks=plot_levels * plot_norm
            )  # Add a colorbar with a nice size
        cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar
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
        if style == None:
            c = ax.pcolormesh(
                data["longitude"],
                data["latitude"],
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
                data["longitude"],
                data["latitude"],
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
        cb.set_label(f"Rainfall ({plot_units})")  # Label the colorbar
        ax.set_title(f"Ensemble standard deviation", size=14)  # This plot's title

        fig.suptitle(
            f"Jurre Brishti cGAN forecast: Valid {valid_time} - {valid_time + timedelta(hours=6)}"
        )  # Overall title
        plt.tight_layout()  # Looks nicer
        plt.show()  # Finally draw the plot


# Plot all ensemble members in the cGAN forecast data at a specified valid time.
# Arguments
#    data                  - An xarray DataSet containing the cGAN rainfall forecasts.
#    valid_time_start_hour - The hour the valid time starts at. Can either be 6, 12, 18 or 0.
def plot_GAN_ensemble(data, valid_time_start_hour, style=None, plot_units="mm/h"):

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

    if plot_units == "mm/h":
        plot_norm = 1
    elif plot_units == "mm/day":
        plot_norm = 24
    elif plot_units == "mm/week":
        plot_norm = 7 * 24
    else:
        print(f"ERROR: Unknown plot units {plot_units}")
        print(f"       Options are 'mm/h', 'mm/day', 'mm/week'.")
        return

    # Use a style other than the default
    if style is not None:
        plot_levels, plot_colours = get_contour_levels(style)

    # Convert the forecast valid time to a datetime.datetime format
    valid_time = datetime64_to_datetime(
        data["fcst_valid_time"][0, valid_time_idx].values
    )

    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15)

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
        ax.add_feature(cfeature.BORDERS, linewidth=1)
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
                data["longitude"],
                data["latitude"],
                data["precipitation"][0, ax_idx, valid_time_idx, :, :] * plot_norm,
                norm=colors.LogNorm(*value_range_precip),
                cmap="YlGnBu",
                transform=ccrs.PlateCarree(),
            )
        else:
            c = ax.contourf(
                data["longitude"],
                data["latitude"],
                data["precipitation"][0, ax_idx, valid_time_idx, :, :] * plot_norm,
                colors=plot_colours,
                levels=plot_levels * plot_norm,
                transform=ccrs.PlateCarree(),
            )
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
