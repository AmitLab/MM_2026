import logging
from pathlib import Path
from typing import Optional
from time_utils import get_now_timestemp_as_string

def set_file_logger(logging_file_path: Path, prefix: Optional[str] = None):
    cur_time_str = get_now_timestemp_as_string()
    prefix = "" if prefix is None else f"{prefix}_"
    new_file_stem = f"{prefix}{logging_file_path.stem}_{cur_time_str}"
    logging_file_path = logging_file_path.with_stem(new_file_stem)
    logging.basicConfig(filename=logging_file_path, level=logging.INFO)
