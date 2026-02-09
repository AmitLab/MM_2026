from datetime import datetime
import pytz
TIME_PATTERN = "%Y_%m_%d__%H_%M_%S"
TIME_PATTERN_LEN = 20

def get_now_timestemp_as_string() -> str:
    return datetime.now(tz=pytz.timezone("Israel")).strftime(TIME_PATTERN)

def get_now_date_as_string() -> str:
    return get_now_timestemp_as_string().split("__")[0]
