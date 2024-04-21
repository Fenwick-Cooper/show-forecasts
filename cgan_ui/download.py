import cfgrib
import sysrsync
import xarray as xr
from os import getenv
from pathlib import Path
from typing import Dict
from argparse import ArgumentParser
from datetime import datetime, timedelta
from dask.diagnostics import ProgressBar
from ecmwf.opendata import Client
from ecmwf.opendata.client import Result
from loguru import logger
from cgan_ui.utils import (
    get_data_store_path,
    get_forecast_data_files,
    get_forecast_data_dates,
    get_ifs_forecast_dates,
    get_cgan_forecast_dates,
)
from cgan_ui.constants import AOI_BBOX
from dsrnngan.test_forecast import gen_cgan_forecast


def standardize_dataset(da: xr.DataArray):
    if "x" in da.dims and "y" in da.dims:
        da = da.rename({"x": "lon", "y": "lat"})
    if "longitude" in da.dims and "latitude" in da.dims:
        da = da.rename({"longitude": "lon", "latitude": "lat"})
    return da


def slice_dataset_by_bbox(ds: xr.Dataset, bbox: list[float]):
    ds = ds.sel(lon=slice(bbox[0], bbox[2]))
    if ds.lat.values[0] < ds.lat.values[-1]:
        ds = ds.sel(lat=slice(bbox[1], bbox[3]))
    else:
        ds = ds.sel(lat=slice(bbox[3], bbox[1]))
    return ds


def read_dataset(file_path: str) -> list[xr.DataArray]:
    try:
        ds = cfgrib.open_datasets(file_path)
    except Exception as err:
        logger.error(f"failed to read {file_path} dataset file with error {err}")
        return None
    if type(ds) is list:
        arrays = []
        for i in range(len(ds)):
            if "number" in ds[i].dims:
                arrays.append(standardize_dataset(ds[i]))
        return xr.combine_by_coords(arrays, compat="override")


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
        result = client.download(
            request=request,
            target=target_file,
        )
    except Exception as err:
        logger.error(
            f"failed to download {model} forecast data for {request['date']} with error {err}"
        )
        return None
    else:
        logger.info(f"downloaded {result.urls[0]} successfully")
        return result


def forecast_files_exist(data_date: datetime.date) -> bool:
    data_dates = get_forecast_data_dates()
    return data_date.strftime("%b %d, %Y") in data_dates


def post_process_ecmwf_grib2_dataset(
    grib2_file_name: str,
    source: str | None = "ecmwf",
    stream: str | None = "enfo",
    re_try_times: int | None = 5,
) -> None:
    logger.info(f"executing post-processing task for {grib2_file_name}")
    store_path = get_data_store_path()
    file_path = f"{store_path}/{source}/{stream}/{grib2_file_name}"
    for _ in range(re_try_times):
        ds = read_dataset(file_path)
        if ds is not None:
            break
    if ds is None:
        logger.error(
            f"failed to read {file_path} after {re_try_times} unsuccessful trials"
        )
    else:
        # save entire raw data file to disk in zarr format
        zarr_file = grib2_file_name.replace("grib2", "zarr")
        write_job = ds.chunk().to_zarr(
            f"{store_path}/raw/{source}/{stream}/{zarr_file}", mode="w", compute=False
        )
        with ProgressBar():
            logger.info(f"writing {grib2_file_name} dataset to {zarr_file}")
            write_job.compute()
        # release memory space occupied by write_job
        write_job = None
        # process area specific masks
        for key in AOI_BBOX.keys():
            logger.info(f"processing {grib2_file_name} mask for {key}")
            sliced = slice_dataset_by_bbox(ds, AOI_BBOX[key])
            out_dir = store_path / "interim" / key / source / stream
            if not out_dir.exists():
                out_dir.mkdir(parents=True)
            sliced.to_netcdf(out_dir / grib2_file_name.replace("grib2", "nc"))
            logger.info(f"saved {grib2_file_name} mask for {key} successfully")
        # remove grib2 file from disk
        # grib2_path: Path = store_path / source / stream / grib2_file_name
        # logger.info(f"deleting {grib2_file_name} from disk")
        # try:
        #     grib2_path.unlink(missing_ok=False)
        # except Exception as err:
        #     logger.error(f"failed to delete {grib2_file_name} with error {err}")
        # # remove associated idx file
        # idx_path: Path = store_path / source / stream / f"{grib2_file_name}.923a8.idx"
        # idx_path.unlink(missing_ok=True)


def post_process_downloaded_ecmwf_forecasts(
    source: str | None = "ecmwf", stream: str | None = "enfo"
) -> None:
    grib2_files = sorted(get_forecast_data_files(source=source, stream=stream))
    logger.info(
        f"starting batch post-processing tasks for {'  <---->  '.join(grib2_files)}"
    )
    for grib2_file in grib2_files:
        post_process_ecmwf_grib2_dataset(
            source=source, stream=stream, grib2_file_name=grib2_file
        )


def syncronize_open_ifs_forecast_data(
    source: str | None = "ecmwf",
    model: str | None = "ifs",
    resolution: str | None = "0p25",
    stream: str | None = "enfo",
    date_str: str | None = None,
    dateback: int | None = 4,
    start_step: int | None = 30,
    final_step: int | None = 54,
    re_try_times: int | None = 10,
    force_download: bool | None = False,
    min_grib2_size: float | None = 4.1 * 1024,
) -> None:
    status_file = Path(getenv("LOGS_DIR", "./")) / "open-ifs.status"

    # check if there is an active data syncronization job
    with open(status_file, "r") as sf:
        status = sf.read()

    # proceed only if there is no active data syncronization job
    if not status:
        # generate down download parameters
        data_dates = get_possible_forecast_dates(data_date=date_str, dateback=dateback)
        steps = get_relevant_forecast_steps(start=start_step, final=final_step)

        # construct data store path
        data_path = get_data_store_path() / source / stream

        # create data directory if it doesn't exist
        if not data_path.exists():
            data_path.mkdir(parents=True, exist_ok=True)

        # create data download client
        client = Client(source=source, model=model, resol=resolution)
        # get latest available forecast date
        latest_fdate = client.latest()

        # set data syncronization status
        with open(status_file, "w") as sf:
            sf.write(1)

        for data_date in data_dates:
            if latest_fdate.date() >= data_date:
                requests = [
                    {
                        "date": data_date,
                        "step": step,
                        "stream": stream,
                    }
                    for step in steps
                ]
                grib2_files = []
                for request in requests:
                    file_name = f"{request['date'].strftime('%Y%m%d')}000000-{request['step']}h-{stream}-ef.grib2"
                    target_file = f"{str(data_path)}/{file_name}"
                    if (
                        not Path(target_file).exists()
                        or Path(target_file).stat().st_size / (1024 * 1024)
                        < min_grib2_size
                        or force_download
                    ):
                        get_url = client._get_urls(
                            request=request, target=target_file, use_index=False
                        )
                        logger.info(
                            f"trying {model} data download with payload {request} on URL {get_url.urls[0]}"
                        )
                        for _ in range(re_try_times):
                            result = try_data_download(
                                client=client,
                                request=request,
                                target_file=target_file,
                                model=model,
                            )
                            if result is not None:
                                logger.info(
                                    f"dataset for {model} forecast, {request['step']}h step, {result.datetime} successfully downloaded"
                                )
                                break
                    else:
                        logger.warning(
                            f"data download job for {request['step']}h {data_date} not executed because the file exist. Pass force_download=True to re-download the files"
                        )
                    grib2_files.append(file_name)
                for grib2_file in grib2_files:
                    post_process_ecmwf_grib2_dataset(
                        source=source, stream=stream, grib2_file_name=grib2_file
                    )
            else:
                logger.warning(
                    f"IFS forecast data for {data_date} is not available. Please try again later!"
                )

        # set data syncronization status
        with open(status_file, "w") as sf:
            sf.write(-1)


def generate_cgan_forecasts():
    ifs_dates = get_ifs_forecast_dates()
    gan_dates = get_cgan_forecast_dates()
    for ifs_date in ifs_dates:
        if ifs_date not in gan_dates:
            # generate forecast for date
            gen_cgan_forecast(
                in_ifs_file=f"IFS_{ifs_date.year}{ifs_date.month:02}{ifs_date.day:02}_00Z.nc"
            )


def syncronize_post_processed_ifs_data(verbose: bool | None = False):
    status_file = Path(getenv("LOGS_DIR", "./")) / "post-processed-ifs.status"

    # check if there is an active data syncronization job
    with open(status_file, "r") as sf:
        status = sf.read()

    if not status:
        src_ssh = getenv("IFS_SERVER", "username@domain.example")
        assert (
            src_ssh != "username@domain.example"
        ), "you must specify IFS data source server address"
        src_dir = getenv("IFS_DIR", "/data/Operational/")
        dest_dir = get_data_store_path() / "IFS"
        logger.info("starting syncronization of IFS forecast data")

        # set data syncronization status
        with open(status_file, "w") as sf:
            sf.write(1)

        sysrsync.run(
            source=src_dir,
            destination=dest_dir,
            source_ssh=src_ssh,
            sync_source_contents=True,
            options=["-a", "-v", "-P"],
            verbose=verbose,
        )

        generate_cgan_forecasts()

        # set data syncronization status
        with open(status_file, "w") as sf:
            sf.write(-1)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-c",
        "--command",
        dest="command",
        type=str,
        default="download",
        help="Command to be executed. Either of download or process.",
    )
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
    dict_args = {key: value for key, value in args.__dict__.items() if key != "command"}
    match (args.command):
        case "download":
            logger.info(
                f"received ecmwf forecast data download task with parameters {dict_args}"
            )
            syncronize_open_ifs_forecast_data(**dict_args)
        case "process":
            logger.info(
                "received ecmwf forecast datasets post-processing task for initial download grib2 files"
            )
            post_process_downloaded_ecmwf_forecasts()
        case _:
            logger.error(f"handler for {args.command} not implemented!")
