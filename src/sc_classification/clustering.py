from typing import Optional, List

from omegaconf import DictConfig
import anndata as ad
import scanpy as sc

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
