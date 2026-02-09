from datetime import date
from pathlib import Path
from typing import Optional


def generate_path_in_output_dir(config, file, add_version=False, add_date_timestamp=False,
                                with_version: Optional[str] = None, with_date_timestamp: Optional[str] = None):
    file_path = Path(config.outputs.output_dir, file)
    if add_version or (with_version is not None):
        if add_version and (with_version is not None):
            raise ValueError("can only get version from either 'add_version' or 'with_version', got from both")
        version = with_version if with_version is not None else config.data_loading.version
        file_path = file_path.with_stem(file_path.stem + f"_data_v_{version}")
    if add_date_timestamp or (with_date_timestamp is not None):
        if add_date_timestamp and (with_date_timestamp is not None):
            raise ValueError(
                "can only get version from either 'add_date_timestamp' or 'with_date_timestamp', got from both")
        date_ts = with_date_timestamp if with_date_timestamp is not None else date.today().isoformat()
        file_path = file_path.with_stem(file_path.stem + f"_ts_{date_ts}")

    return file_path
