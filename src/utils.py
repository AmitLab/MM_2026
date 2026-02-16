import pandas as pd
import numpy as np
import os
import subprocess
from typing import Optional, Union, Dict, Any
import matplotlib.axes
import scanpy as sc
from anndata import AnnData
from matplotlib.colors import Colormap, TwoSlopeNorm
from scanpy.plotting._utils import savefig_or_show
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns
import psutil
import warnings
from sklearn.preprocessing import StandardScaler

from gtfparse import read_gtf
from rnanorm import TPM
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def calculate_features_ms(exp_tpm, genesets_ms_dict, way = 'median'):
    y = StandardScaler().fit_transform(exp_tpm)
    df_y = pd.DataFrame(y, columns=exp_tpm.columns, index=exp_tpm.index)

    l = []
    feats = ["CD1", "CD2", "PR", "MS", "MF", "HP", "LB"]
    for feat in feats:
        if way == 'mean':
            _ = df_y.T.reindex(
                genesets_ms_dict[feat + "_UP"]
            ).dropna(axis=0).mean(axis=0) - df_y.T.reindex(
                genesets_ms_dict[feat + "_DN"]
            ).dropna(
                axis=0
            ).mean(
                axis=0
            )
        elif way == 'median':
            _ = df_y.T.reindex(
                genesets_ms_dict[feat + "_UP"]
            ).dropna(axis=0).median(axis=0) - df_y.T.reindex(
                genesets_ms_dict[feat + "_DN"]
            ).dropna(
                axis=0
            ).median(
                axis=0
            )
        else:
            warnings.warn("Choose either 'mean' or 'median' for the 'way' parameter.")
            return None
        l.append(_)
    Features = pd.concat(l, axis=1)
    Features.columns = feats
    return Features

# batch aware scaling for MS

def calculate_features_ms_percentile(exp_tpm_log, genesets_ms_dict, way='median'):
    df_percentile = exp_tpm_log.rank(axis=0, pct=True)
    
    l = []
    feats = ["CD1", "CD2", "PR", "MS", "MF", "HP", "LB"]
    for feat in feats:
        up_genes = df_percentile.T.reindex(genesets_ms_dict[feat + "_UP"]).dropna(axis=0)
        dn_genes = df_percentile.T.reindex(genesets_ms_dict[feat + "_DN"]).dropna(axis=0)
        
        if way == 'mean':
            up_score = up_genes.mean(axis=0)
            dn_score = dn_genes.mean(axis=0)
        elif way == 'median':
            up_score = up_genes.median(axis=0)
            dn_score = dn_genes.median(axis=0)
        else:
            import warnings
            warnings.warn("Choose either 'mean' or 'median' for the 'way' parameter.")
            return None
        
        score = up_score - dn_score
        l.append(score)
    
    Features = pd.concat(l, axis=1)
    Features.columns = feats
    return Features

def intersect_series_with_index(index: pd.Index, data: pd.Series) -> pd.Series:
    aligned_series = data.loc[index.intersection(data.index)]
    return aligned_series

def copy_data(source_path, destination_path, filename):
    subprocess.run(['cp', os.path.join(source_path, filename), 
                    os.path.join(destination_path, filename)], 
                   check=True)
    
def prepare_xrep_adata(adata: AnnData, use_rep: str = "X_cnv", cmap: str = 'bwr', **kwargs) -> (AnnData, Dict[str, Any]):

    chr_pos_dict = dict(sorted(adata.uns['cnv']["chr_pos"].items(), key=lambda x: x[1]))
    chr_pos = list(chr_pos_dict.values())

    tmp_adata = AnnData(X=adata.obsm[use_rep], obs=adata.obs, uns=adata.uns)

    tmp_data = tmp_adata.X.data if issparse(tmp_adata.X) else tmp_adata.X
    vmin = kwargs.pop("vmin", None)
    vmax = kwargs.pop("vmax", None)
    if vmin is None:
        vmin = np.nanmin(tmp_data)
    if vmax is None:
        vmax = np.nanmax(tmp_data)
    kwargs["norm"] = TwoSlopeNorm(0, vmin=vmin, vmax=vmax)

    # add chromosome annotations
    var_group_positions = list(zip(chr_pos, chr_pos[1:] + [tmp_adata.shape[1]]))
    tmp_adata.uns['var_group_positions'] = var_group_positions
    tmp_adata.uns['var_group_labels'] = list(chr_pos_dict.keys())

    return tmp_adata

def intersect_index(index: pd.Index, data: pd.Series) -> pd.Series:
    aligned_data = data.loc[index.intersection(data.index)]
    return aligned_data

def item_series(item, indexed=None):
    if indexed is not None:
        if hasattr(indexed, 'index'):
            return pd.Series([item] * len(indexed), index=indexed.index)
        elif type(indexed) is int and indexed > 0:
            return pd.Series([item] * indexed, index=np.arange(indexed))
    return pd.Series()

def subsample_adata_per_category(adata, n_cells_per_population, category='Populations', random_state=0):

    np.random.seed(random_state)
    subsampled_indices = []

    for population in adata.obs[category].unique():
        population_indices = np.where(adata.obs[category] == population)[0]
        
        if len(population_indices) > n_cells_per_population:
            subsampled_population_indices = np.random.choice(population_indices, n_cells_per_population, replace=False)
        else:
            subsampled_population_indices = population_indices
        
        subsampled_indices.extend(subsampled_population_indices)

    subsampled_adata = adata[subsampled_indices].copy()
    
    return subsampled_adata

def check_memory():
    process = psutil.Process()
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss / (1024 ** 2):.2f} MB")

def load_intervals(data_path='/home/labs/amit/annaku/repos/MM_2024_AK/data/',
                   filename='gene_intervals_hg38P.ens90.txt'):
    intervals = pd.read_csv(os.path.join(data_path, 
                            filename), 
                            index_col = 'gene_name', 
                            sep = '\t').rename(columns={'chrom': 'chromosome'})
    return intervals.groupby(intervals.index).agg({
        'chromosome': 'first',
        'start': 'first',
        'end': 'last',
        'strand': 'first'
    })

def calculate_cdf(s: pd.Series):
    return s.sort_values().cumsum() / s.sum()

def median_scale(data, mask=None, axis=0):
    if mask is not None:
        exclude_idx = mask[~mask].index
        ref = data.reindex(data.index.intersection(exclude_idx))
    else:
        ref = data

    center = ref.median(axis=axis)

    if isinstance(data, pd.Series):
        spread = mad(ref.dropna())
        if spread == 0:
            spread = 1.0
        result = (data - center) / spread
    else:
        broadcast_axis = (axis + 1) % 2
        spread = ref.apply(lambda col: mad(col.dropna()), axis=axis)
        spread = spread.replace(0, 1.0)
        result = data.sub(center, axis=broadcast_axis).div(spread, axis=broadcast_axis)

    return result

def weights_multipl_run(exp, coef_df):
    # coed_df - dataframe, where indexes are genes, columns are features, values are weights
    # exp - log transformed tpm exp, indexes are samples, columns are genes
    gene_names_common = coef_df.index.intersection(exp.columns)
    score = exp[gene_names_common] @ coef_df.loc[gene_names_common]
    return score

def intersect_df(dataframes=()):
    shared_indices = set(dataframes[0].index)
    for df_index in range(1, len(dataframes)):
        shared_indices = shared_indices.intersection(dataframes[df_index].index)
    return [dataframes[df_index].loc[list(shared_indices)] for df_index in range(len(dataframes))]

def format_pvalue(p, digits=3, use_stars=False, prefix='p-value'):
    significant_p = 10**-digits
    if use_stars:
        pvalue_str = star_pvalue(p, lev3=10**-digits)
    else:
        if p < significant_p:
            if len(prefix):
                prefix += ' < '
            pvalue_str = f'{prefix}{significant_p}'
        else:
            if len(prefix):
                prefix += ' = '
            if p < 0.00001:
                pvalue_str = f'{prefix}{round_to_1(p):.0e}'
            else:
                pvalue_str = f'{prefix}{round_to_1(p)}'
    return pvalue_str

def assign_quantiles(data, quantiles=[0.5]):
    quantiles = sorted(set([0] + quantiles + [1]))
    quantile_values = data.quantile(quantiles)
    labels = [f'{lower}q<x<={upper}q' for lower, upper in zip(quantiles[:-1], quantiles[1:])]
    return pd.cut(data, bins=quantile_values, labels=labels, include_lowest=True)

def categorize_high_low(data, quantiles=[0.5]):
    return assign_quantiles(data, quantiles).map(
        {f'0q<x<{str(quantiles[0])}q': 'Low', f'{str(quantiles[0])}q<x<1q': 'High'}
    )

def generate_color_palette(series, cmap = plt.cm.gist_rainbow, min_v=0., max_v=1.):

    unique_values = series.sort_values().unique()
    n_colors = len(unique_values)
    colors = plt.cm.get_cmap(cmap)(np.linspace(min_v, max_v, n_colors))
    hsv_colors = matplotlib.colors.rgb_to_hsv(colors[:, :3])
    rgb_colors = matplotlib.colors.hsv_to_rgb(hsv_colors)
    hex_colors = np.apply_along_axis(matplotlib.colors.to_hex, 1, rgb_colors)
    return dict(zip(unique_values, hex_colors))

def load_gmt(path):
    gene_signatures = {}
    with open(path, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            name = parts[0]
            genes = parts[2:]
            gene_signatures[name] = genes
    return gene_signatures

def calculate_tmm(counts_df):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['counts'] = counts_df
        
        ro.r('''
        library(edgeR)
        dge <- DGEList(counts = counts)
        dge <- calcNormFactors(dge, method = "TMM", logratioTrim=.3, sumTrim=0.05)
        norm_counts <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)
        ''')
        
        norm_counts = ro.globalenv['norm_counts']
    
    return pd.DataFrame(norm_counts, columns=counts_df.columns, index=counts_df.index)

## calculate_tmm is the same as following code in the notebook:
# %%R -i counts -o norm_counts -o gene_names -o sample_names
# library(edgeR)
# dge <- DGEList(counts = counts)
# dge <- calcNormFactors(dge, method = "TMM",
# logratioTrim=.3, sumTrim=0.05, # default values
# )
# norm_counts <- cpm(dge, log = FALSE, normalized.lib.sizes = TRUE)
# gene_names <- rownames(norm_counts)
# sample_names <- colnames(norm_counts)

def calculate_tpm_tmm_layers(pdata, target_genes, gtf_path):
    pdata_tm = pdata[:, [gene for gene in target_genes if gene in pdata.var_names]].copy()
    exp = pdata_tm.to_df()
    # TPM
    dataset = Dataset(exp=exp.copy(), gtf_path=gtf_path)
    tpm = TPM(dataset.gtf_path).set_output(transform="pandas")
    exp_tpm = tpm.fit_transform(dataset.exp)
    exp_tpm = exp_tpm.rename(columns=dataset.gene_mapping_reverse)
    # TMM
    exp_renamed = exp.rename(columns=dataset.gene_mapping_reverse)
    counts = exp_renamed.T
    norm_counts_df = calculate_tmm(counts)
    norm_counts_df = np.log2(norm_counts_df + 1)
    expr_my_tmm = norm_counts_df.T
    
    exp_tpm_aligned = exp_tpm.loc[pdata_tm.obs_names, pdata_tm.var_names]
    pdata_tm.layers['TPM'] = exp_tpm_aligned.values
    pdata_tm.layers['log1p_TPM'] = np.log1p(exp_tpm_aligned.values)
    
    expr_my_tmm_aligned = expr_my_tmm.loc[pdata_tm.obs_names, pdata_tm.var_names]
    pdata_tm.layers['TMM'] = expr_my_tmm_aligned.values
    pdata_tm.layers['log1p_TMM'] = np.log1p(expr_my_tmm_aligned.values)
    
    return pdata_tm

class Dataset:
    """Helper class for TPM calculation with gene ID mapping"""
    def __init__(self, exp, gtf_path):
        self.exp = exp
        self.gtf_path = gtf_path
        self.gene_mapping_reverse = None
        
        from gtfparse import read_gtf
        gtf_df = read_gtf(self.gtf_path)
        gtf_df = pd.DataFrame(gtf_df, columns=gtf_df.columns)
        gene_info = gtf_df[gtf_df['feature'] == 'gene'][['gene_name', 'gene_id']]
        gene_mapping = dict(zip(gene_info['gene_name'], gene_info['gene_id']))
        self.gene_mapping_reverse = dict(zip(gene_info['gene_id'], gene_info['gene_name']))
        exp.rename(columns=gene_mapping, inplace=True)