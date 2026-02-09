import anndata as ad
import pandas as pd


def add_number_of_patients_in_neighborhood(adata: ad.AnnData, patient_col: str,
                                           number_of_patient_col="number_of_diffrent_patients_in_nighborhood",
                                           return_as_series=False, distances_obsp_key="distances"):
    coo_distance = adata.obsp[distances_obsp_key].tocoo(copy=True)
    cids = adata.obs_names
    cell_patients_count = pd.DataFrame({
        "cells": cids[coo_distance.col],
        "neigbors_cell_patient": cids[coo_distance.row].map(adata.obs[patient_col])
    }).groupby('cells').neigbors_cell_patient.nunique()
    adata.obs[number_of_patient_col] = cell_patients_count
    if return_as_series:
        return cell_patients_count


def count_number_of_annotation_in_neighborhood(adata, annotation_col, annotation_val,
                                               annotation_count_col=None, distances_obsp_key="distances"):
    annotation_count_col = annotation_count_col if annotation_count_col is not None else f"count_of_{annotation_val}_in_neighborhood"
    coo_distance = adata.obsp[distances_obsp_key].tocoo(copy=True)
    cids = adata.obs_names
    coo_annotations = pd.DataFrame({
        "cell": cids[coo_distance.col],
        "neighbors_annotation": cids[coo_distance.row].map(adata.obs[annotation_col])
    })
    coo_annotations["neighbors_annotation"] = coo_annotations["neighbors_annotation"] == annotation_val

    adata.obs[annotation_count_col] = coo_annotations.groupby('cell').neighbors_annotation.sum()
