import logging
import os
from abc import ABC, abstractmethod
from functools import partial
from pathlib import Path
from typing import Callable, List, Optional

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from tqdm.contrib.concurrent import process_map

from data_loading.utils import load_dataframe_from_file


class AnnDataLoader(ABC):

    @abstractmethod
    def load_data_to_anndata(self) -> ad.AnnData:
        pass


DEBUG_MODE = False
DEBUG_N_BATCHES = 10


class MultiMethodFromPlatesDataLoader(AnnDataLoader):
    def __init__(self, sc_data_dirs: List[Path], plates_data_paths: List[Path], plate_id_col_name: str,
                 plates_data_transform_functions: Optional[List[Callable[[pd.DataFrame], pd.DataFrame]]] = None):
        if len(sc_data_dirs) != len(plates_data_paths):
            raise ValueError(f"got sc_data_dirs and plates_data_paths with different length, both need to much")
        self.loaders = []
        for sc_data_dir, plates_data_path in zip(sc_data_dirs, plates_data_paths):
            self.loaders.append(FromPlatesDataLoader(sc_data_dir=sc_data_dir,
                                                     plates_data_path=plates_data_path,
                                                     plate_id_col_name=plate_id_col_name,
                                                     plates_data_transform_functions=plates_data_transform_functions))

    def load_data_to_anndata(self) -> ad.AnnData:
        ads = [loader.load_data_to_anndata() for loader in self.loaders]
        return ad.concat(adatas=ads)


class FromPlatesDataLoader(AnnDataLoader):
    def __init__(self, sc_data_dir: Path, plates_data_path: Path, plate_id_col_name: str,
                 plates_data_transform_functions: Optional[List[Callable[[pd.DataFrame], pd.DataFrame]]] = None):
        self.sc_data_dir = sc_data_dir
        self.plates_data_path = plates_data_path
        self.plate_id_col_name = plate_id_col_name
        self.plates_data_transform_functions = plates_data_transform_functions if plates_data_transform_functions is not None else []

    def _get_single_plate(self, row_tpl, col_names):
        row = row_tpl[1]
        logging.info(f"Reading , plate id - {row[self.plate_id_col_name]}")
        cur_data = sc.read_text(Path(self.sc_data_dir, row[self.plate_id_col_name] + ".txt"))
        cur_data = cur_data.T
        logging.debug("inserting metadata to the anndata")
        for col_name in col_names:
            cur_data.obs[col_name] = row[col_name]
        logging.debug(f"converting data from batch id - {row[self.plate_id_col_name]} to sparse matrix")
        cur_data.X = csr_matrix(cur_data.X)
        return cur_data

    def load_data_to_anndata(self) -> ad.AnnData:
        plates_data_df = load_dataframe_from_file(self.plates_data_path)
        plates_data_df = plates_data_df.dropna(axis=0, subset=[self.plate_id_col_name])
        if DEBUG_MODE:
            self.plates_data_transform_functions.append(lambda df: df.head(DEBUG_N_BATCHES))
        for transform_func in self.plates_data_transform_functions:
            logging.info(f"applying the function: {transform_func.__name__} on the meta data")
            plates_data_df = transform_func(plates_data_df)

        # Read all plates into anndata and merge them
        col_names = plates_data_df.columns
        if DEBUG_MODE:
            adatas = [self._get_single_plate(batch_tpl, col_names) for batch_tpl in plates_data_df.iterrows()]
        else:
            adatas = process_map(partial(self._get_single_plate, col_names=col_names),
                                 list(plates_data_df.iterrows()), max_workers=os.cpu_count() // 2,
                                 desc="loading relevant plates", chunksize=5,
                                 unit="plate")
        logging.info("merging to single adata")
        adata = ad.concat(adatas, merge="same")

        adata.obs = adata.obs.astype(str)

        return adata
