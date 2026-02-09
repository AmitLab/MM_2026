import logging
from typing import List

import anndata as ad
import scanpy as sc

from pre_processing.black_list import get_mask_that_match_any_of_patterns


def pp_drop_cells(adata, min_num_genes_per_cell: int, min_num_counts_per_cell: int,
                  mitochondrial_patterns: List[str], max_pct_mt_genes_pre_cell: float):
    logging.info("dropping bad cells")
    n_cells_before = len(adata.obs)

    sc.pp.filter_cells(adata, min_genes=min_num_genes_per_cell)
    n_cells_after_min_gene_count = len(adata.obs)
    logging.info(f"dropped cells by min genes , dropped: {n_cells_before - n_cells_after_min_gene_count}")

    sc.pp.filter_cells(adata, min_counts=min_num_counts_per_cell)

    n_cells_after_filter = len(adata.obs)
    logging.info(f"dropped cells by min counts, dropped: {n_cells_after_min_gene_count - n_cells_after_filter}")

    adata.var['mt'] = get_mask_that_match_any_of_patterns(adata.var_names, mitochondrial_patterns)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    mt_subset_index = adata.obs.index[(adata.obs.pct_counts_mt < max_pct_mt_genes_pre_cell)]
    adata._inplace_subset_obs(mt_subset_index)
    n_cells_after_mt = len(adata.obs)
    logging.info(f"dropped cells by mitochondrial percentage, dropped: {n_cells_after_filter - n_cells_after_mt}")

    n_cells_dropped = n_cells_before - n_cells_after_mt
    logging.info(f"dropped over all, before: {n_cells_before}, after: {n_cells_after_mt}, dropped: {n_cells_dropped}")

# def pp_drop_cells_other(adata, min_num_genes_per_cell: int, min_num_counts_per_cell: int,
#                   mitochondrial_patterns: List[str], max_pct_mt_genes_pre_cell: float):
#     logging.info("dropping bad cells")
#     n_cells_before = len(adata.obs)

#     data.var['mt'] = get_mask_that_match_any_of_patterns(adata.var_names, mitochondrial_patterns)
#     sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#     sc.pp.filter_cells(adata, min_genes=min_num_genes_per_cell)
#     n_cells_after_min_gene_count = len(adata.obs)
#     logging.info(f"dropped cells by min genes , dropped: {n_cells_before - n_cells_after_min_gene_count}")

#     sc.pp.filter_cells(adata, min_counts=min_num_counts_per_cell)

#     n_cells_after_filter = len(adata.obs)
#     logging.info(f"dropped cells by min counts, dropped: {n_cells_after_min_gene_count - n_cells_after_filter}")

#     # adata.var['mt'] = get_mask_that_match_any_of_patterns(adata.var_names, mitochondrial_patterns)
#     # sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#     mt_subset_index = adata.obs.index[(adata.obs.pct_counts_mt < max_pct_mt_genes_pre_cell)]
#     adata._inplace_subset_obs(mt_subset_index)
#     n_cells_after_mt = len(adata.obs)
#     logging.info(f"dropped cells by mitochondrial percentage, dropped: {n_cells_after_filter - n_cells_after_mt}")

#     n_cells_dropped = n_cells_before - n_cells_after_mt
#     logging.info(f"dropped over all, before: {n_cells_before}, after: {n_cells_after_mt}, dropped: {n_cells_dropped}")


def pp_drop_genes(adata: ad.AnnData, min_num_cells_per_gene):
    logging.info("dropping genes by min number of cells ")
    n_genes_before = len(adata.var)

    sc.pp.filter_genes(adata, min_cells=min_num_cells_per_gene)

    n_genes_after = len(adata.var)
    n_genes_dropped = n_genes_before - n_genes_after
    logging.info(f"dropped genes, before: {n_genes_before}, after: {n_genes_after}, dropped: {n_genes_dropped}")
