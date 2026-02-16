from typing import Optional, List, Iterable

import anndata as ad
import pandas as pd
import scanpy as sc
from omegaconf import DictConfig

from pre_processing.black_list import drop_blacklist_genes_by_pattern


def normalize_and_choose_genes(adata: ad.AnnData, conf: DictConfig, target_sum=1e4,
                               genes_to_keep: Optional[List[str]] = None, **highly_variable_genes_kwargs) -> ad.AnnData:
    adata_for_clustering = adata.copy()

    sc.pp.normalize_total(adata_for_clustering, target_sum=target_sum)
    sc.pp.log1p(adata_for_clustering)
    adata_for_clustering.raw = adata_for_clustering  # freeze the state in `.raw`

    if genes_to_keep is None:
        sc.pp.highly_variable_genes(adata_for_clustering, **highly_variable_genes_kwargs, subset=True)
    else:
        adata_for_clustering = adata_for_clustering[:, genes_to_keep]

    adata_for_clustering = drop_blacklist_genes_by_pattern(adata_for_clustering,
                                                           conf.sc_classification.gene_selection_patterns_blacklist).copy()

    return adata_for_clustering


def shuang_genes_to_keep(genes_names: Iterable[str], flavor: str) -> Optional[List[str]]:
    if flavor == "None":
        return None
    elif flavor == 'MARS_SPID_common':
        common_mars_spid_genes_shuang = pd.read_csv('/home/projects/amit/annaku/data/mm_2023/feats/common_genes.csv')
        allowed_genes = list(common_mars_spid_genes_shuang["0"])
    elif flavor == 'MARS_SPID_combined':
        combined_mars_spid_genes_shuang = pd.read_csv('/home/projects/amit/annaku/data/mm_2023/feats/combined_genes.csv')
        allowed_genes = list(combined_mars_spid_genes_shuang["0"])
    else:
        raise ValueError("value in flavor not supported,"
                         " supported are: 'None', 'MARS_SPID_common' or 'MARS_SPID_combined'")
    genes_to_keep = [gene for gene in genes_names if gene in allowed_genes]
    return genes_to_keep
