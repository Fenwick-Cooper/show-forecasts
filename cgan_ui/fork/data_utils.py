# Utility functions and data used internally for processing forecasts
import json
from typing import Dict
from pathlib import Path
import numpy as np
from os import getenv
from cgan_ui.constants import COUNTRY_NAMES
from datetime import datetime, timedelta, timezone
import shapefile
import cartopy.io.shapereader as shpreader


def get_locations_data() -> list[Dict[str, str]]:
    data_file = f"{getenv('APP_DIR', '.')}/shapefiles/locations.json"
    with open(data_file, "r") as jf:
        return json.loads(jf.read())


def get_shape_boundary(
    shape_name: str | None = COUNTRY_NAMES[0],
) -> shpreader.BasicReader:
    # Get the shapefile
    try:
        return shpreader.Reader(f"{getenv('APP_DIR', '.')}/shapefiles/{shape_name}.shp")
    except Exception:
        return shpreader.Reader(
            f"{getenv('APP_DIR', '.')}/shapefiles/{COUNTRY_NAMES[0]}.shp"
        )


# Returns the normalisation used to plot
# Arguments
#   plot_units='mm/h' - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
# Returns
#   The normalisation to apply when plotting or 1 if the units are not specified correctly.
#   The string plot_units corresponding to the normalisation.
def get_plot_normalisation(plot_units):
    if plot_units == "mm/h":
        plot_norm = 1
    elif plot_units == "mm/6h":
        plot_norm = 6
    elif plot_units == "mm/day":
        plot_norm = 24
    elif plot_units == "mm/week":
        plot_norm = 7 * 24
    else:
        print(f"ERROR: Unknown plot units '{plot_units}'.")
        print(f"       Options are 'mm/h', 'mm/6h', 'mm/day', 'mm/week'.")
        print(f"       Selecting 'mm/h'.")
        plot_norm = 1
        plot_units = "mm/h"
    return plot_norm, plot_units


# Returns the bounding box of the region that we want to plot
# Arguments
#   region='ICPAC' - can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                    'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#    border_size   - Area around the region in degrees to include in the plot
def get_region_extent(
    shape_name: str | None = COUNTRY_NAMES[0], border_size: float | None = 0.5
):
    try:
        sf = shapefile.Reader(f"{getenv('APP_DIR', '.')}/shapefiles/{shape_name}.shp")
    except Exception:
        sf = shapefile.Reader(
            f"{getenv('APP_DIR', '.')}/shapefiles/{COUNTRY_NAMES[0]}.shp"
        )
    # find boundary index
    if (
        shape_name != COUNTRY_NAMES[0]
        and not Path(f"{getenv('APP_DIR', '.')}/shapefiles/{shape_name}.shp").exists()
    ):
        shape_index = [
            index
            for index in range(len(sf.records()))
            if sf.record(index).as_dict()["name"] == shape_name
        ]
        bbox = None if not len(shape_index) else sf.shape(shape_index[0]).bbox
    else:
        bbox = sf.bbox

    if bbox is None:
        print(
            f"ERROR: Boundary name {shape_name} is not recognized. Should be one of {', '.join(COUNTRY_NAMES)}"
        )
        return None
    # The boundary of the region in the form of a cartopy extent
    return [
        bbox[0] - border_size,
        bbox[2] + border_size,
        bbox[1] - border_size,
        bbox[3] + border_size,
    ]


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
        plot_levels = np.array([0, 50, 100, 200, 1000]) / (
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

        plot_levels = np.array([0, 1, 10, 30, 50, 100, 200, 1000]) / (
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
            np.array([0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 1000])
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
        plot_levels = np.array([0, 50, 100, 200, 1000]) / (
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
            np.array([0, 1, 5, 10, 20, 30, 50, 75, 100, 1000]) / 24
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


# Convert numpy.datetime64 to datetime
# Arguments
#    datetime64 - A single date stored in a datetime64 object.
# Returns
#    A date as a datetime.datetime object.
def datetime64_to_datetime(datetime64):
    unix_epoch = np.datetime64(0, "s")
    one_second = np.timedelta64(1, "s")
    seconds_since_epoch = (datetime64 - unix_epoch) / one_second
    time_zone = timezone(timedelta(0), name="UTC")
    return datetime.fromtimestamp(seconds_since_epoch, time_zone)
    # return datetime.utcfromtimestamp(seconds_since_epoch)


# Prints the named locations available for making histograms
# Arguments
#   country=None - Optionally restrict the country to one of 'Kenya', 'South Sudan', 'Rwanda',
#                  'Burundi', 'Djibouti', 'Eritrea', 'Ethiopia', 'Sudan', 'Somalia',
#                  'Tanzania', 'Uganda'
def print_locations(country=None):
    with open("shapefiles/locations.json") as jf:
        locations = json.loads(jf.read())

    for i in range(len(locations)):
        if (country == None) or (country == locations[i]["country"]):
            print(
                f"{locations[i]['name']}, {locations[i]['country']}, ({locations[i]['latitude']}N, {locations[i]['longitude']}E)"
            )
