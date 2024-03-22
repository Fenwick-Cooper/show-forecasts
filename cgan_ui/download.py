from ecmwf.opendata import Client
from datetime import datetime, timedelta
from loguru import logger
from pathlib import Path
import os


def get_possible_forecast_dates(dateback: int | None = 4) -> list[datetime.date]:
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


def download_ifs_forecast_data(
    source: str | None = "ecmwf",
    model: str | None = "ifs",
    resolution: str | None = "0p25",
    stream: str | None = "enfo",
    dateback: int | None = 4,
    start_step: int | None = 30,
    final_step: int | None = 54,
) -> None:
    # generate down download parameters
    dates = get_possible_forecast_dates(dateback)
    steps = get_relevant_forecast_steps(start=start_step, final=final_step)
    client = Client(source=source, model=model, resol=resolution)
    data_store = os.getenv("APP_STORE_DIR", str(Path("./store")))
    requests = [
        {
            "date": date,
            "step": step,
            "stream": stream,
        }
        for date in dates
        for step in steps
    ]

    # create data directory if it doesn't exist
    data_path = Path(f"{data_store}/{source}/{stream}")
    if not data_path.exists():
        data_path.mkdir(parents=True)

    for request in requests:
        target_file = f"{data_path}/{request['date'].strftime('%Y%m%d')}000000-{request['step']}h-{stream}-ef.grib2"
        get_url = client._get_urls(request=request, target=target_file, use_index=False)
        logger.info(
            f"trying {model} data download with payload {request} on URL {get_url.urls[0]}"
        )
        try:
            result = client.download(
                request=request,
                target=target_file,
            )
        except Exception as err:
            logger.error(
                f"failed to download {model} forecast data for {request['date']} with error {err}"
            )
        else:
            logger.info(
                f"downloaded {model} forecast data for {result.datetime} successfully"
            )


if __name__ == "__main__":
    download_ifs_forecast_data()
