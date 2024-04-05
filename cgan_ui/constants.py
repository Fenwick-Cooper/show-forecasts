# Required information about each variable we want to look at
DATA_PARAMS = {
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
LEAD_START_HOUR = 30
LEAD_END_HOUR = 54
