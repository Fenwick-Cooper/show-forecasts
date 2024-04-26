# Required information about each variable we want to look at
DATA_PARAMS = {
    "tp": {
        "name": "Total precipitation",
        "units": "mm/day",
        "normalisation": 0.001,  # Convert from m to mm/day
        "offset": 0,
        "accumulated": True,
    },
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
LEAD_START_HOUR = 30
LEAD_END_HOUR = 54

# area of interest bounding box
AOI_BBOX = {
    "EA": [21, -11.75, 51, 24],
    "KEN": [33.6, -5.84, 43.6, 6.22],
    "ETH": [31.58, 3.12, 47.65, 16.5],
}

# visualization color schemes
COLOR_SCHEMES = ["ICPAC", "ICPAC_heavy", "KMD", "EMI", "EMI_heavy", "Default"]

# boundary layers
REGION_SHAPES = ["ICPAC", "Kenya", "Ethiopia"]

TP_PLOT_UNITS = ["mm/week", "mm/day", "mm/6h", "mm/h"]
