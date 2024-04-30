from enum import Enum


class DataSourceIdentifier(Enum):
    ecmwf = "open_ifs"
    gbmc = "gbmc_ifs"
    cgan = "cgan_forecast"
