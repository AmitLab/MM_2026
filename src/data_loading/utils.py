from functools import partial
from pathlib import Path
from typing import List, Optional

import anndata as ad
import pandas as pd

_supported_suffixes_to_reader = {
    ".csv": pd.read_csv,
    ".xlsx": pd.read_excel,
    ".txt": partial(pd.read_csv, delimiter="\t"),
    ".tsv": partial(pd.read_csv, delimiter="\t")
}


def load_dataframe_from_file(plates_data_path: Path) -> pd.DataFrame:
    path_suffix = plates_data_path.suffix
    if path_suffix in _supported_suffixes_to_reader:
        plates_data_df = _supported_suffixes_to_reader[path_suffix](plates_data_path)
        return plates_data_df
    else:
        raise ValueError(f"Not supperted meta_data_file_type:"
                         f" got {path_suffix}, supported are {list(_supported_suffixes_to_reader)}")


def merge_labels_to_adata(adata: ad.AnnData, labels_df: pd.DataFrame, col_in_adata_to_merge_by: str,
                          cols_in_labels_df_to_merge_by: str, cols_to_validate_not_empty: Optional[List[str]],
                          labels_col_names_to_merge: Optional[List[str]] = None,
                          merge_suffixes=("_x", "_y")) -> ad.AnnData:
    new_adata = adata.copy()
    if cols_in_labels_df_to_merge_by == "index":
        labels_df = labels_df.reset_index()
    if col_in_adata_to_merge_by == "index":
        new_adata.obs.index = new_adata.obs.index.map(str)
    else:
        new_adata.obs = new_adata.obs.astype({col_in_adata_to_merge_by: "str"})

    if labels_col_names_to_merge is None:
        labels_col_names_to_merge = list(labels_df.columns)
    if cols_in_labels_df_to_merge_by not in labels_col_names_to_merge:
        labels_col_names_to_merge.append(cols_in_labels_df_to_merge_by)

    new_adata.obs = pd.merge(
        left=new_adata.obs.reset_index(),
        right=labels_df[labels_col_names_to_merge], how="left", validate="m:1", suffixes=merge_suffixes,
        left_on=col_in_adata_to_merge_by, right_on=cols_in_labels_df_to_merge_by).set_index("index")
    if cols_to_validate_not_empty is not None:
        new_adata = new_adata[~ new_adata.obs[cols_to_validate_not_empty].isna().any(axis=1), :].copy()
    return new_adata


def get_sample_level_data(adata_path: Path, clinical_data_path: Path, architypes_df: pd.DataFrame,
                          sample_creation_kwargs: dict, clinical_loading_kwargs: dict):
    adata = ad.read_h5ad(adata_path)
    adata.obs['Hospital.Code'] = adata.obs['Hospital.Code'].apply(
        lambda h_code: f"cart_p{h_code.replace('cart', '')}" if h_code in ('cart13', 'cart21') else h_code)
    samples_df = extract_samples_metadata(adata, **sample_creation_kwargs)


def extract_samples_metadata(adata: ad.AnnData, metadata_cols, split_by_method=True, split_by_sample=True,
                             generate_architype_id=True, generate_hl_architype_id=True, generate_sample_level_ids=False,
                             code_lower_case=True):
    groupby_cols = ['Hospital.Code']
    if split_by_method:
        groupby_cols.append('Method')

    if split_by_sample:
        groupby_cols.append('Biopsy.Sequence')
    else:
        metadata_cols.append('Biopsy.Sequence')

    metadata_cols = [col for col in metadata_cols if col not in groupby_cols]
    all_samples = adata.obs.groupby(groupby_cols)[metadata_cols].agg(lambda x: pd.Series.mode(x).iloc[0]).dropna()
    all_samples.columns.name = None
    all_samples = all_samples.reset_index()
    if code_lower_case:
        all_samples['Hospital.Code'] = all_samples['Hospital.Code'].str.lower()

    if generate_architype_id:
        if not split_by_method:
            raise ValueError(
                "cannot create architype id with seperation of samples by method, use split_by_method=True")
        all_samples['PID'] = 'z.' + all_samples['Method'].astype(str) + '_malignant_' + all_samples[
            'Hospital.Code'].astype(str)
        all_samples['PID'] = all_samples['PID'].str.lower()
        if generate_sample_level_ids:
            all_samples['SID'] = 'z.' + all_samples['Method'].astype(str) + '_malignant_' + \
                                 all_samples['Hospital.Code'].astype(str) + '_' + all_samples['Biopsy.Sequence'].astype(
                str)
            all_samples['SID'] = all_samples['SID'].str.lower()
    if generate_hl_architype_id:
        if not split_by_method:
            raise ValueError(
                "cannot create architype id with seperation of samples by method, use split_by_method=True")
        all_samples['PID_Healthy_Like'] = 'z.' + all_samples['Method'].astype(str) + '_healthy_like_' + all_samples[
            'Hospital.Code'].astype(str)
        all_samples['PID_Healthy_Like'] = all_samples['PID_Healthy_Like'].str.lower()
        if generate_sample_level_ids:
            all_samples['SID_Healthy_Like'] = 'z.' + all_samples['Method'].astype(str) + '_healthy_like_' + \
                                              all_samples['Hospital.Code'].astype(str) + '_' + all_samples[
                                                  'Biopsy.Sequence'].astype(str)
            all_samples['SID_Healthy_Like'] = all_samples['SID_Healthy_Like'].str.lower()

    all_samples['Biopsy.Sequence'] = all_samples['Biopsy.Sequence'].astype(int)

    return all_samples


def get_updated_disease_col(metadata_df: pd.DataFrame, disease_col: str, hospital_disease_col: str,
                            update_non_naive_NDMM_by_treatment_history: bool, remove_PRMM: bool,
                            treatment_names: Optional[List[str]] = None,
                            non_naive_NDMM_value="non_naive_NDMM") -> pd.Series:
    new_disease_col = metadata_df.apply(
        lambda row: row[disease_col] if pd.isna(row[hospital_disease_col]) else row[hospital_disease_col], axis=1)

    if update_non_naive_NDMM_by_treatment_history:
        if treatment_names is None:
            treatment_names = ["Bortezomib", "Ixazomib", "Carfilzomib", "Lenalidomide", "Thalidomide", "Pomalidomide",
                               "Cyclophosphamide", "Chemotherapy", "Venetoclax", "Dexamethasone", "Prednisone",
                               "Daratumumab", "Elotuzumab", "Belantamab", "Talquetamab", "Teclistamab", "Cevostamab",
                               "Selinexor", "Auto-SCT", "CART"]
        non_naive_NDMM_mask = (new_disease_col == "NDMM") & (metadata_df[treatment_names].any(axis=1))
        new_disease_col[non_naive_NDMM_mask] = non_naive_NDMM_value

    if remove_PRMM:
        new_disease_col = new_disease_col.apply(lambda x: "RRMM" if x == "PRMM" else x)

    return new_disease_col
