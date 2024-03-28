import os
from pathlib import Path
from typing import Dict
from argparse import ArgumentParser
from datetime import datetime, timedelta
from ecmwf.opendata import Client
from ecmwf.opendata.client import Result
from loguru import logger
from cgan_ui.show_forecasts import var_info, get_forecast_data_dates


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


def try_data_download(
    client: Client,
    request: Dict[str, str],
    target_file: str,
    model: str | None = "ifs",
) -> Result | None:
    try:
        result = client.retrieve(
            request=request,
            target=target_file,
        )
    except Exception as err:
        logger.error(
            f"failed to download {model} forecast data for {request['date']} with error {err}"
        )
        return None
    else:
        return result


def forecast_files_exist(data_date: datetime.date) -> bool:
    data_dates = get_forecast_data_dates()
    return data_date.strftime("%b %d, %Y") in data_dates


def download_ifs_forecast_data(
    source: str | None = "ecmwf",
    model: str | None = "ifs",
    resolution: str | None = "0p25",
    stream: str | None = "enfo",
    date_str: str | None = None,
    dateback: int | None = 4,
    start_step: int | None = 30,
    final_step: int | None = 54,
    re_try_times: int | None = 5,
    force_download: bool | None = False,
) -> None:
    # generate down download parameters
    data_dates = get_possible_forecast_dates(data_date=date_str, dateback=dateback)
    steps = get_relevant_forecast_steps(start=start_step, final=final_step)

    # construct data store path
    data_store = os.getenv("DATA_STORE_DIR", str(Path("./store")))
    # create data directory if it doesn't exist
    data_path = Path(f"{data_store}/{source}/{stream}")
    if not data_path.exists():
        data_path.mkdir(parents=True, exist_ok=True)

    # create data download client
    client = Client(source=source, model=model, resol=resolution)
    # get latest available forecast date
    latest_fdate = client.latest()

    for data_date in data_dates:
        if latest_fdate.date() >= data_date:
            requests = [
                {
                    "date": data_date,
                    "step": step,
                    "param": list(var_info.keys()),
                    "stream": stream,
                }
                for step in steps
            ]
            for request in requests:
                target_file = f"{data_path}/{request['date'].strftime('%Y%m%d')}000000-{request['step']}h-{stream}-ef.grib2"
                if (
                    not Path(target_file).exists()
                    or Path(target_file).stat().st_size / (1024 * 1024) < 80
                    or force_download
                ):
                    get_url = client._get_urls(
                        request=request, target=target_file, use_index=False
                    )
                    logger.info(
                        f"trying {model} data download with payload {request} on URL {get_url.urls[0]}"
                    )
                    for i in range(re_try_times):
                        result = try_data_download(
                            client=client,
                            request=request,
                            target_file=target_file,
                            model=model,
                        )
                        if result is not None:
                            logger.info(
                                f"downloaded {model} forecast data for {request['step']}h {result.datetime} successfully"
                            )
                            break
                else:
                    logger.warning(
                        f"data download job for {request['step']}h {data_date} not executed because the file exist. Pass force_download=True to re-download the files"
                    )
        else:
            logger.warning(
                f"IFS forecast data for {data_date} is not available. Please try again later!"
            )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--date",
        dest="date_str",
        type=str,
        default=None,
        help="Download ECMWF forecast data for given date",
    )
    parser.add_argument(
        "-p",
        "--period",
        dest="dateback",
        type=int,
        default=4,
        help="Forecasts data time range in days since current date",
    )
    parser.add_argument(
        "-s",
        "--start",
        dest="start_step",
        type=int,
        default=30,
        help="Forecast data start time step",
    )
    parser.add_argument(
        "-f",
        "--final",
        dest="final_step",
        type=int,
        default=54,
        help="Forecast data final time step",
    )
    args = parser.parse_args()
    download_ifs_forecast_data(**args.__dict__)
