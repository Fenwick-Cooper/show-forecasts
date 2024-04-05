import schedule
import time
from datetime import datetime
from cgan_ui.download import download_ifs_forecast_data
from loguru import logger

for hour in range(8, 23, 1):
    schedule.every().day.at(f"{str(hour).rjust(2, '0')}:00", "Africa/Nairobi").do(
        download_ifs_forecast_data, date_str=datetime.now().strftime("%Y-%m-%d")
    )

for job in schedule.get_jobs():
    logger.info(f"scheduled ecmwf data download task {job}")

while True:
    all_jobs = schedule.get_jobs()
    schedule.run_pending()
    time.sleep(30 * 60)
