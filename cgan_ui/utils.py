import os
from datetime import datetime
from pathlib import Path


def get_data_store_path() -> Path:
    return Path(os.getenv("DATA_STORE_DIR", str(Path("./store")))).absolute()


def get_forecast_data_files(
    mask: str | None = None, source: str | None = "ecmwf", stream: str | None = "enfo"
) -> list[str]:
    store_path = get_data_store_path()
    data_dir = (
        store_path / source / stream
        if mask is None
        else store_path / "interim" / mask / source / stream
    )
    return [str(dfile).split("/")[-1] for dfile in data_dir.iterdir()]


def get_forecast_data_dates(
    mask: str | None = None, source: str | None = "ecmwf", stream: str | None = "enfo"
) -> list[str]:
    data_files = get_forecast_data_files(source=source, stream=stream, mask=mask)
    # extract forecast initialization dates
    data_dates = sorted(sorted(set([dfile.split("-")[0] for dfile in data_files])))
    return [
        datetime.strptime(dt, "%Y%m%d%H%M%S").strftime("%b %d, %Y") for dt in data_dates
    ]
