import sys, os, time, schedule
from pathlib import Path
from datetime import datetime
from cgan_ui.download import (
    syncronize_open_ifs_forecast_data,
    syncronize_post_processed_ifs_data,
)
from loguru import logger

logger_opts = dict(
    enqueue=True,
    backtrace=True,
    diagnose=True,
    format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}",
)

config = {
    "handlers": [
        {"sink": sys.stdout, "colorize": True, **logger_opts},
        {
            "sink": Path(os.getenv("LOGS_DIR", "./")) / "cgan-jobs.log",
            "serialize": True,
            **logger_opts,
        },
    ],
    "extra": {"app": "cgan-jobs"},
}
logger.configure(**config)


for hour in range(12, 00, 1):
    schedule.every().day.at(f"{str(hour).rjust(2, '0')}:00", "Africa/Nairobi").do(
        syncronize_open_ifs_forecast_data, date_str=datetime.now().strftime("%Y-%m-%d")
    )
    schedule.every().day.at(f"{str(hour).rjust(2, '0')}:00", "Africa/Nairobi").do(
        syncronize_post_processed_ifs_data
    )

for job in schedule.get_jobs():
    logger.info(f"scheduled ecmwf data download task {job}")

while True:
    all_jobs = schedule.get_jobs()
    schedule.run_pending()
    time.sleep(30 * 60)
