from functools import partial
from pathlib import Path

import pandas as pd

_supported_suffixes_to_reader = {
    ".csv": pd.read_csv,
    ".xlsx": pd.read_excel,
    ".txt": partial(pd.read_csv, delimiter="\t")
}


def load_plates_data_from_file(plates_data_path: Path) -> pd.DataFrame:
    path_suffix = plates_data_path.suffix
    if path_suffix in _supported_suffixes_to_reader:
        plates_data_df = _supported_suffixes_to_reader[path_suffix](plates_data_path)
        return plates_data_df
    else:
        raise ValueError(f"Not supperted meta_data_file_type:"
                         f" got {path_suffix}, supported are {list(_supported_suffixes_to_reader)}")
