import sysrsync, subprocess, cfgrib
import xarray as xr
from os import getenv
from pathlib import Path
from typing import Dict
from argparse import ArgumentParser
from datetime import datetime
from ecmwf.opendata import Client
from ecmwf.opendata.client import Result
from loguru import logger
from cgan_ui.utils import (
    get_data_store_path,
    get_forecast_data_files,
    get_ifs_forecast_dates,
    get_cgan_forecast_dates,
    get_possible_forecast_dates,
    get_relevant_forecast_steps,
    get_data_sycn_status,
    set_data_sycn_status,
)
from cgan_ui.constants import AOI_BBOX


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


def post_process_ecmwf_grib2_dataset(
    grib2_file_name: str,
    source: str | None = "ecmwf",
    stream: str | None = "enfo",
    re_try_times: int | None = 5,
    force_process: bool | None = False,
) -> None:
    logger.info(f"executing post-processing task for {grib2_file_name}")
    store_path = get_data_store_path()
    ea_nc_file = (
        store_path
        / "interim"
        / "EA"
        / source
        / stream
        / grib2_file_name.replace("grib2", "nc")
    )
    if not ea_nc_file.exists() or force_process:
        logger.info(f"post-processing ECMWF IFS forecast data file {grib2_file_name}")
        grib2_dir = store_path / source / stream
        file_path = f"{grib2_dir}/{grib2_file_name}"
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
            # zarr_file = grib2_file_name.replace("grib2", "zarr")
            # write_job = ds.chunk().to_zarr(
            #     f"{store_path}/raw/{source}/{stream}/{zarr_file}",
            #     mode="w",
            #     compute=False,
            # )
            # with ProgressBar():
            #     logger.info(f"writing {grib2_file_name} dataset to {zarr_file}")
            #     write_job.compute()
            # # release memory space occupied by write_job
            # write_job = None
            # process area specific masks
            for key in AOI_BBOX.keys():
                logger.info(f"processing {grib2_file_name} mask for {key}")
                sliced = slice_dataset_by_bbox(ds, AOI_BBOX[key])
                out_dir = store_path / "interim" / key / source / stream
                if not out_dir.exists():
                    out_dir.mkdir(parents=True)
                nc_file = out_dir / grib2_file_name.replace("grib2", "nc")
                if nc_file.exists():
                    nc_file.unlink()

                sliced.to_netcdf(out_dir / grib2_file_name.replace("grib2", "nc"))
                logger.info(f"saved {grib2_file_name} mask for {key} successfully")
            # remove grib2 file from disk
            logger.info(f"archiving {grib2_file_name}")
            archive_dir = store_path / "archive" / source / stream
            if not archive_dir.exists():
                archive_dir.mkdir(parents=True)

            try:
                Path(file_path).replace(target=archive_dir / f"{grib2_file_name}")
            except Exception as err:
                logger.error(f"failed to archive {grib2_file_name} with error {err}")
            # remove idx files from the disk
            idx_files = [
                idxf for idxf in grib2_dir.iterdir() if idxf.name.endswith(".idx")
            ]
            for idx_file in idx_files:
                try:
                    idx_file.unlink()
                except Exception as err:
                    logger.error(
                        f"failed to delete grib2 index file {idx_file.name} with error {err}"
                    )


def post_process_downloaded_ecmwf_forecasts(
    source: str | None = "ecmwf", stream: str | None = "enfo"
) -> None:
    grib2_files = sorted(get_forecast_data_files(source=source, stream=stream))
    logger.info(
        f"starting batch post-processing tasks for {'  <---->  '.join(grib2_files)}"
    )
    for grib2_file in grib2_files:
        post_process_ecmwf_grib2_dataset(
            source=source, stream=stream, grib2_file_name=grib2_file, force_process=True
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
    min_nc_size: float | None = 360 * 1024,
) -> None:
    logger.info(
        f"recived IFS open forecast data syncronization job at {datetime.now().strftime('%Y-%m-%d %H:%M')}"
    )
    # proceed only if there is no active data syncronization job
    if not get_data_sycn_status():
        logger.info(
            f"starting IFS open forecast data syncronization at {datetime.now().strftime('%Y-%m-%d %H:%M')}"
        )
        # generate down download parameters
        data_dates = get_possible_forecast_dates(data_date=date_str, dateback=dateback)
        steps = get_relevant_forecast_steps(start=start_step, final=final_step)

        # construct data store path
        data_path = get_data_store_path() / source / stream
        mask_path = get_data_store_path() / "interim" / "EA" / source / stream

        # create data directory if it doesn't exist
        if not data_path.exists():
            data_path.mkdir(parents=True, exist_ok=True)

        # create data download client
        client = Client(source=source, model=model, resol=resolution)
        # get latest available forecast date
        latest_fdate = client.latest()

        # set data syncronization status
        set_data_sycn_status(source="ecmwf", status=1)

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
                    mask_file = mask_path / file_name.replace(".grib2", ".nc")
                    target_file = data_path / file_name
                    target_size = (
                        0
                        if not target_file.exists()
                        else target_file.stat().st_size / (1024 * 1024)
                    )
                    mask_size = (
                        0
                        if not mask_file.exists()
                        else mask_file.stat().st_size / (1024 * 1024)
                    )
                    if (
                        not (target_file.exists() or mask_file.exists())
                        or target_size < min_grib2_size
                        or mask_size < min_nc_size
                        or force_download
                    ):
                        get_url = client._get_urls(
                            request=request, target=str(target_file), use_index=False
                        )
                        logger.info(
                            f"trying {model} data download with payload {request} on URL {get_url.urls[0]}"
                        )
                        for _ in range(re_try_times):
                            result = try_data_download(
                                client=client,
                                request=request,
                                target_file=str(target_file),
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
        set_data_sycn_status(source="ecmwf", status=0)


def generate_cgan_forecasts():
    ifs_dates = get_ifs_forecast_dates()
    gan_dates = get_cgan_forecast_dates()
    for ifs_date in ifs_dates:
        if ifs_date not in gan_dates:
            # generate forecast for date
            in_ifs_file = (
                f"IFS_{ifs_date.year}{ifs_date.month:02}{ifs_date.day:02}_00Z.nc"
            )
            subprocess.call(
                shell=True,
                cwd=f'{getenv("WORK_HOME","/opt/mycgan")}/ensemble-cgan/dsrnngan',
                args=f"python test_forecast.py -f {in_ifs_file}",
            )


def syncronize_post_processed_ifs_data(verbose: bool | None = False):
    logger.info(
        f"recived post processed IFS forecast data syncronization job at {datetime.now().strftime('%Y-%m-%d %H:%M')}"
    )
    if not get_data_sycn_status(source="cgan"):
        logger.info(
            f"starting post processed IFS forecast data syncronization at {datetime.now().strftime('%Y-%m-%d %H:%M')}"
        )
        ifs_host = getenv("IFS_SERVER_HOST", "domain.example")
        ifs_user = getenv("IFS_SERVER_USER", "username")
        src_ssh = f"{ifs_user}@{ifs_host}"
        assert (
            src_ssh != "username@domain.example"
        ), "you must specify IFS data source server address"
        src_dir = getenv("IFS_DIR", "/data/Operational")
        dest_dir = get_data_store_path() / "IFS"
        logger.info("starting syncronization of IFS forecast data")

        # set data syncronization status
        set_data_sycn_status(source="cgan", status=1)

        sysrsync.run(
            source=str(src_dir),
            destination=f"{str(dest_dir)}/",
            source_ssh=src_ssh,
            private_key=getenv("IFS_PRIVATE_KEY", "/srv/ssl/private.key"),
            sync_source_contents=True,
            options=["-a", "-v", "-P"],
            verbose=verbose,
        )

        generate_cgan_forecasts()

        # set data syncronization status
        set_data_sycn_status(source="cgan", status=0)


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
        case "cgan":
            syncronize_post_processed_ifs_data()
        case _:
            logger.error(f"handler for {args.command} not implemented!")
