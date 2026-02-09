import pandas as pd
import scanpy as sc
import logging
import infercnvpy as cnv
import anndata as ad
from src.utils import *
import numpy as np
from pathlib import Path

from typing import Optional, Union, Dict, Any

import matplotlib.axes
from anndata import AnnData
from matplotlib.colors import Colormap, TwoSlopeNorm
from scanpy.plotting._utils import savefig_or_show
from scipy.sparse import issparse


def summarize_group(group, cells_number_to_sum=50):
    num_cells = len(group)
    num_minibulk = num_cells // cells_number_to_sum  
    summarized_data = []
    
    for i in range(num_minibulk):
        start_idx = i * cells_number_to_sum
        end_idx = (i + 1) * cells_number_to_sum
        minibulk_indices = np.arange(start_idx, end_idx)
        
        summed_counts = group.X[minibulk_indices].sum(axis=0).A1  # Ensure 1D array
        
        minibulk_obs = group.obs.iloc[start_idx].copy()
        minibulk_obs['n_cells'] = len(minibulk_indices)
        
        summarized_data.append((summed_counts, minibulk_obs))
    
    return summarized_data

def standart_pp(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return adata

def create_minibulk_and_filter(adata, cells_number_to_sum=50, cat_to_group = ['Hospital.Code', 'leiden_healthy'], to_filter = True, total_counts_thr2 = 155000,
                               preproc_standart = True):
    grouped = adata.obs.groupby(cat_to_group)
    summarized_data = grouped.apply(lambda x: summarize_group(adata[x.index], cells_number_to_sum = cells_number_to_sum))
    summarized_data = [item for sublist in summarized_data for item in sublist]

    minibulk_counts = np.array([x[0] for x in summarized_data])
    minibulk_obs = pd.DataFrame([x[1] for x in summarized_data])
    adata_minibulk = ad.AnnData(X=minibulk_counts, obs=minibulk_obs, var=adata.var)
    
    sc.pp.calculate_qc_metrics(adata_minibulk, percent_top=None, log1p=True, inplace=True)
    
    if to_filter == True:
        mask = adata_minibulk.obs['total_counts'] < total_counts_thr2 #to ensure 
        adata_minibulk = adata_minibulk[mask]

    adata_minibulk.raw = adata_minibulk
    
    if preproc_standart == True:
        adata_minibulk = standart_pp(adata_minibulk)

    return adata_minibulk

def load_intervals(data_path = '/home/projects/amit/annaku/repos/MM_2024_AK/data/'):
    intervals = pd.read_csv(os.path.join(data_path, 
                            'gene_intervals_hg38P.ens90.txt'), 
                            index_col = 'gene_name', 
                            sep = '\t').rename(columns={'chrom': 'chromosome'})
    return intervals.groupby(intervals.index).agg({
        'chromosome': 'first',
        'start': 'first',
        'end': 'last',
        'strand': 'first'
    })

def process_minibulk_cnv(adata,
                         batch = 'Method',
                         cat_to_group = ['Sample.Code', 'leiden_healthy'],
                         reference_key = 'leiden_healthy',
                         reference_cat=['Healthy'],
                         to_filter = True,
                         preproc_standart = True,
                         window_size=400):

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    minibulk_list = []

    for method in adata.obs[batch].unique():
        adata_minibulk_method = create_minibulk_and_filter(
            adata[adata.obs[batch] == method], 
            cells_number_to_sum=50, 
            cat_to_group=cat_to_group, 
            to_filter=to_filter, 
            total_counts_thr2=155000,
            preproc_standart=preproc_standart
        )
        logger.info(f"Created minibulk for method: {method}")
        
        adata_minibulk_method.var = adata_minibulk_method.var.merge(load_intervals(), left_index=True, right_index=True)
        logger.info(f"Added coordinates to minibulk for method: {method}")
        
        cnv.tl.infercnv(
            adata_minibulk_method, 
            reference_key=reference_key, 
            reference_cat=reference_cat,
            window_size=window_size,
            step=1
        )
        logger.info(f"Ran infercnv for minibulk of method: {method}")
        
        minibulk_list.append(adata_minibulk_method)

    adata_minibulk_combined = ad.concat(minibulk_list, join='outer', uns_merge='unique')
    logger.info("Combined minibulks into a single AnnData object")

    cnv.tl.pca(adata_minibulk_combined)
    cnv.pp.neighbors(adata_minibulk_combined)
    cnv.tl.leiden(adata_minibulk_combined)
    cnv.tl.umap(adata_minibulk_combined)
    # cnv.tl.cnv_score(adata_minibulk_combined)
    logger.info("Processed combined minibulk data")

    return adata_minibulk_combined

# adata_minibulk_combined = process_minibulk_cnv(adata_cnv)

# adata_minibulk_combined.write_h5ad(os.path.join(data_path, f'adata_minibulk_cnv_v_{version}_ws_{window_size}_only_pc_annotated_filtered.h5ad'))

def prepare_xrep_adata(adata: AnnData, use_rep: str = "X_cnv", cmap: str = 'bwr', **kwargs) -> (AnnData, Dict[str, Any]):


    chr_pos_dict = dict(sorted(adata.uns['cnv']["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())

    q_band_pos_dict = dict(sorted(adata.uns['cnv']["q_band_pos"].items(), key=lambda x: x[1]))
    q_band_pos = list(q_band_pos_dict.values())

    var = pd.DataFrame(index=[f"{i}" for i in range(adata.obsm[use_rep].shape[1])])
    tmp_adata = AnnData(X=adata.obsm[use_rep], obs=adata.obs, var=var, uns=adata.uns)

    tmp_data = tmp_adata.X.data if issparse(tmp_adata.X) else tmp_adata.X
    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    if vmin is None:
        vmin = np.nanmin(tmp_data)
    if vmax is None:
        vmax = np.nanmax(tmp_data)
    kwargs["norm"] = TwoSlopeNorm(0, vmin=vmin, vmax=vmax) # original norm, I will use another

    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [tmp_adata.shape[1]]))

    var_group_positions_band = list(zip(q_band_pos, q_band_pos[1:] + [tmp_adata.shape[1]]))

    kwargs['cmap'] = cmap
    tmp_adata.uns['var_group_positions_band'] = var_group_positions_band
    tmp_adata.uns['var_group_labels_band'] = list(q_band_pos_dict.keys())
    tmp_adata.uns['q_band_pos'] = q_band_pos

    tmp_adata.uns['var_group_positions'] = var_group_positions
    tmp_adata.uns['var_group_labels'] = list(chr_pos_dict.keys())
    tmp_adata.uns['chr_pos'] = chr_pos
    #tmp_adata.uns['q_band_pos'] = q_band_pos

    return tmp_adata


def average_X_per_category(adata, category1='Sample.Code', category2=None):
    if category2:
        combined_category = adata.obs[category1].astype(str) + "_" + adata.obs[category2].astype(str)
    else:
        combined_category = adata.obs[category1]
    
    unique_categories = combined_category.unique()  # preserve order
    
    averaged_X = np.zeros((len(unique_categories), adata.X.shape[1]))
    obs_first = []
    
    for i, cat in enumerate(unique_categories):
        mask = combined_category == cat
        averaged_X[i, :] = np.mean(adata.X[mask, :], axis=0)
        obs_first.append(adata.obs[mask].iloc[0])
    
    obs_first = pd.DataFrame(obs_first, columns=adata.obs.columns)
    
    averaged_adata = AnnData(
        X=averaged_X, 
        obs=obs_first, 
        var=adata.var, 
        uns=adata.uns, 
        obsm=adata.obsm, 
        obsp=adata.obsp
    )
    
    return averaged_adata

def sliding_window_cnv_arms(adata, window_size=10, window_step=5, min_genes=20):

    arm_data = []
    arm_names = []
    arm_metrics = [] 
    
    unique_chr_arms = adata.var['chr_arm'].dropna().unique()
    
    for chr_arm in sorted(unique_chr_arms):
        arm_genes = adata.var.index[adata.var['chr_arm'] == chr_arm]
        
        if len(arm_genes) >= min_genes:
            arm_data_subset = adata[:, arm_genes].X
            
            windows = []
            for i in range(0, len(arm_genes) - window_size + 1, window_step):
                window_data = arm_data_subset[:, i:i+window_size]
                window_mean = np.mean(window_data, axis=1)
                windows.append(window_mean)
            
            windows_array = np.array(windows).T  

            arm_mean = np.mean(windows_array, axis=1)
            arm_std = np.std(windows_array, axis=1)
            arm_max = np.max(windows_array, axis=1)
            arm_min = np.min(windows_array, axis=1)
            
            arm_data.append(arm_mean)
            arm_names.append(chr_arm)
            arm_metrics.append({
                'n_genes': len(arm_genes),
                'n_windows': len(windows),
                'std': arm_std,
                'max': arm_max,
                'min': arm_min
            })
    
    X = np.column_stack(arm_data)
    var_df = pd.DataFrame(arm_metrics, index=arm_names)
    
    adata_arms = ad.AnnData(
        X=X,
        obs=adata.obs,
        var=var_df
    )
    
    return adata_arms

def sliding_window_cnv_region(adata, region_col = 'chr_arm', window_size=10, window_step=5, min_genes=20):

    # regions defined in adata.var column

    region_data = []
    region_names = []
    region_metrics = [] 
    
    unique_regions = adata.var[region_col].dropna().unique()
    
    for region in sorted(unique_regions):
        region_genes = adata.var.index[adata.var[region_col] == region]
        
        if len(region_genes) >= min_genes:
            region_data_subset = adata[:, region_genes].X
            
            windows = []
            for i in range(0, len(region_genes) - window_size + 1, window_step):
                window_data = region_data_subset[:, i:i+window_size]
                window_mean = np.mean(window_data, axis=1)
                windows.append(window_mean)
            
            windows_array = np.array(windows).T 

            region_mean = np.mean(windows_array, axis=1)
            region_std = np.std(windows_array, axis=1)
            region_max = np.max(windows_array, axis=1)
            region_min = np.min(windows_array, axis=1)
            
            region_data.append(region_mean)
            region_names.append(region)
            region_metrics.append({
                'n_genes': len(region_genes),
                'n_windows': len(windows),
                'std': region_std,
                'max': region_max,
                'min': region_min
            })
    
    X = np.column_stack(region_data)
    var_df = pd.DataFrame(region_metrics, index=region_names)
    
    adata_region = ad.AnnData(
        X=X,
        obs=adata.obs,
        var=var_df
    )
    
    return adata_region

def sort_genes_by_position(
    adata: ad.AnnData,
    chr_col: str = 'chr',
    start_col: str = 'start'
):
    var_df = adata.var.copy()
    var_df[chr_col] = var_df[chr_col].astype(str).str.replace('chr', '', regex=False)
    
    def safe_to_int(x):
        try:
            return int(x)
        except ValueError:
            return x  
    
    var_df[chr_col] = var_df[chr_col].apply(safe_to_int)
    var_df = var_df.sort_values(by=[chr_col, start_col], ascending=True)
    return adata[:, var_df.index].copy()


def calculate_genome_wide_cnv_burden(
    adata: ad.AnnData,
    window_size: int = 100,
    window_step: int = 25,
    baseline: float = 1.0,
    sort_genes: bool = True,
    chr_col: str = 'chr',
    start_col: str = 'start'
):
    if sort_genes:
        adata = sort_genes_by_position(adata, chr_col, start_col)

    X = adata.X 
    n_patients, n_genes = X.shape
    burdens = np.zeros(n_patients, dtype=float)
    stds = np.zeros(n_patients, dtype=float)
    
    for start_idx in range(0, n_genes - window_size + 1, window_step):
        window_data = X[:, start_idx : start_idx + window_size]
        window_mean = window_data.mean(axis=1)
        window_std = window_data.std(axis=1)
        burdens += np.abs(window_mean - baseline)
        stds += window_std
    adata.obs['cnv_burden'] = burdens
    adata.obs['cnv_burden_std'] = stds
    
    return adata