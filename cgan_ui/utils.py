import os, re
from datetime import datetime, timedelta
from pathlib import Path

from cgan_ui.constants import LEAD_START_HOUR, LEAD_END_HOUR, DATA_PARAMS


def get_data_store_path() -> Path:
    return Path(os.getenv("DATA_STORE_DIR", str(Path("./store")))).absolute()


def get_possible_forecast_dates(
    data_date: str | None = None, dateback: int | None = 4
) -> list[datetime.date]:
    if data_date is not None:
        return [datetime.strptime(data_date, "%Y-%m-%d").date()]
    now = datetime.now()
    dates = [now.date()]
    for i in range(1, dateback + 1):
        new_date = now - timedelta(days=i)
        dates.append(new_date.date())
    return dates


def get_relevant_forecast_steps(
    start: int | None = 30, final: int | None = 54, step: int | None = 3
) -> list[int]:
    return [step for step in range(start, final + 1, step)]


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
    mask: str | None = None,
    source: str | None = "ecmwf",
    stream: str | None = "enfo",
    check_steps: bool | None = False,
) -> list[str]:
    data_files = get_forecast_data_files(source=source, stream=stream, mask=mask)
    # extract forecast initialization dates
    match source:
        case "ecmwf":
            # be sure all IFS time steps were synced
            if check_steps:
                r = re.compile(
                    r"([0-9]{4})([0-9]{2})([0-9]{2})000000-([0-9]{2})h-enfo-ef.nc"
                )
                synced_dates = {}
                for fname in data_files:
                    if bool(r.match(fname)):
                        parts = r.split(fname)
                        fdate = f"{parts[1]}{parts[2]}{parts[3]}000000"
                        if fdate in synced_dates.keys():
                            synced_dates[fdate].append(int(parts[4]))
                        else:
                            synced_dates[fdate] = [int(parts[4])]
                print("synced_dates --> ", synced_dates)
                data_dates = list(
                    sorted(
                        [
                            data_date
                            for data_date in synced_dates.keys()
                            if list(sorted(synced_dates[data_date]))
                            == get_relevant_forecast_steps()
                        ]
                    )
                )
                print("data_dates --> ", data_dates)
            else:
                data_dates = sorted(
                    sorted(set([dfile.split("-")[0] for dfile in data_files]))
                )
            return list(
                reversed(
                    [
                        datetime.strptime(dt, "%Y%m%d%H%M%S").strftime("%b %d, %Y")
                        for dt in data_dates
                    ]
                )
            )
        case "cgan":
            return list(
                reversed(
                    [
                        datetime.strptime(dfile, "GAN_%Y%m%d.nc").strftime("%b %d, %Y")
                        for dfile in sorted(data_files)
                    ]
                )
            )
        case _:
            return []


def get_ifs_forecast_dates():
    ifs_dir = get_data_store_path() / "IFS"
    if not ifs_dir.exists():
        ifs_dir.mkdir(parents=True)
    ifs_files = [fpath.name for fpath in ifs_dir.iterdir()]
    return [datetime.strptime(fname, "IFS_%Y%m%d_00Z.nc") for fname in ifs_files]


def get_cgan_forecast_dates():
    cgan_dir = get_data_store_path() / "GAN_forecasts"
    if not cgan_dir.exists():
        cgan_dir.mkdir(parents=True)
    gan_files = [fpath.name for fpath in cgan_dir.iterdir()]
    return [datetime.strptime(fname, "GAN_%Y%m%d.nc") for fname in gan_files]


def set_data_sycn_status(source: str | None = "ecmwf", status: int | None = 1):
    file_name = "post-processed-ifs.log" if source == "cgan" else "ecmwf-open-ifs.log"
    status_file = Path(os.getenv("LOGS_DIR", "./")) / file_name
    # initialize with not-active status
    with open(status_file, "w") as sf:
        sf.write(str(status))


def get_data_sycn_status(source: str | None = "ecmwf") -> int:
    file_name = "post-processed-ifs.log" if source == "cgan" else "ecmwf-open-ifs.log"
    status_file = Path(os.getenv("LOGS_DIR", "./")) / file_name

    if not status_file.exists():
        # initialize with not-active status
        with open(status_file, "w") as sf:
            sf.write("0")

    # check if there is an active data syncronization job
    with open(status_file, "r") as sf:
        return int(sf.read())


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