import os
from datetime import datetime, timedelta
from pathlib import Path

from cgan_ui.constants import LEAD_START_HOUR, LEAD_END_HOUR, DATA_PARAMS


def get_data_store_path() -> Path:
    return Path(os.getenv("DATA_STORE_DIR", str(Path("./store")))).absolute()


def get_forecast_data_files(
    mask: str | None = None, source: str | None = "ecmwf", stream: str | None = "enfo"
) -> list[str]:
    store_path = get_data_store_path()
    match source:
        case "ecmwf":
            data_dir = (
                store_path / source / stream
                if mask is None
                else store_path / "interim" / mask / source / stream
            )
            return [str(dfile).split("/")[-1] for dfile in data_dir.iterdir()]
        case "cgan":
            data_dir = store_path / "GAN_forecasts"
            return [str(dfile).split("/")[-1] for dfile in data_dir.iterdir()]
        case _:
            return []


def get_forecast_data_dates(
    mask: str | None = None, source: str | None = "ecmwf", stream: str | None = "enfo"
) -> list[str]:
    data_files = get_forecast_data_files(source=source, stream=stream, mask=mask)
    # extract forecast initialization dates
    match source:
        case "ecmwf":
            data_dates = sorted(
                sorted(set([dfile.split("-")[0] for dfile in data_files]))
            )
            return [
                datetime.strptime(dt, "%Y%m%d%H%M%S").strftime("%b %d, %Y")
                for dt in data_dates
            ]
        case "cgan":
            return [
                datetime.strptime(dfile, "GAN_%Y%m%d.nc").strftime("%b %d, %Y")
                for dfile in sorted(data_files)
            ]
        case _:
            return []


def get_ifs_forecast_dates():
    ifs_dir = get_data_store_path() / "IFS"
    ifs_files = [fpath.split("/")[-1] for fpath in ifs_dir.iterdir()]
    return [datetime.strptime(fname, "IFS_%Y%m%d_00Z.nc") for fname in ifs_files]


def get_cgan_forecast_dates():
    cgan_dir = get_data_store_path() / "GAN_forecasts"
    gan_files = [fpath.split("/")[-1] for fpath in cgan_dir.iterdir()]
    return [datetime.strptime(fname, "GAN_%Y%m%d.nc") for fname in gan_files]


# Some info to be clear with dates and times
# Arguments
#    forecast_init_date - A datetime.datetime corresponding to when the forecast was initialised.
def print_forecast_info(forecast_init_date):
    start_date = forecast_init_date + timedelta(hours=LEAD_START_HOUR)
    end_date = forecast_init_date + timedelta(hours=LEAD_END_HOUR)
    print(f"Forecast average: {start_date} - {end_date}")
    print(f"Forecast initialisation: {forecast_init_date.date()} 00:00:00")
    print()


# Lists the variables available for plotting
# Arguments
#    print_variables=True - print a list of the variables (True or False)
# Returns
#    A list of the variable names that can be used for plotting.
def get_possible_variables(print_variables=True):

    keys = list(DATA_PARAMS.keys())
    if print_variables:
        print("Available variables to plot are the following:")
        for i in range(len(keys)):
            print(
                f"{keys[i]:<5} - {DATA_PARAMS[keys[i]]['name']} ({DATA_PARAMS[keys[i]]['units']})"
            )
        print()

    return list(DATA_PARAMS.keys())
