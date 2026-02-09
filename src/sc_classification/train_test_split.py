from typing import Optional

import anndata as ad
from sklearn.model_selection import train_test_split


def split_adata(adata: ad.AnnData, split_by_obs_column: Optional[str] = None, train_col_name: str = "train", seed=None) ->\
        (ad.AnnData, ad.AnnData):
    if split_by_obs_column is None:
        train_obs_index, test_obs_index = train_test_split(adata.obs.index, random_state=seed)
        return adata[train_obs_index], adata[test_obs_index]
    else:
        specific_columns_values = set(adata.obs[split_by_obs_column].values)
        train_values, test_values = train_test_split(list(specific_columns_values), random_state=seed)
        adata.obs[train_col_name] = adata.obs[split_by_obs_column].apply(lambda x: x in train_values).astype(bool)
        train_adata = adata[adata.obs.train].copy()
        test_adata = adata[~ adata.obs.train].copy()
        return train_adata, test_adata
