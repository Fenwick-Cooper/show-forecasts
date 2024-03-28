import schedule
import time
from datetime import datetime
from cgan_ui.download import download_ifs_forecast_data

for hour in range(3, 23, 1):
    schedule.every().day.at(f"{str(hour).rjust(2, '0')}:00", "Europe/London").do(
        download_ifs_forecast_data, date_str=datetime.now().strftime("%Y-%m-%d")
    )

while True:
    all_jobs = schedule.get_jobs()
    schedule.run_pending()
    time.sleep(30 * 60)
