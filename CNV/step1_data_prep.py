import os
import sys
import pandas as pd
import numpy as np
import matplotlib
from pathlib import Path

repo_dir = '/home/projects/amit/annaku/repos/Blueprint'
sys.path.append(repo_dir)

from omegaconf import OmegaConf
import scanpy as sc
import anndata as ad

from src.utils import *
from src.cnv_utils import *

version = '20250306'

project_root = '/home/projects/amit/annaku/repos/Blueprint'
conf_path = os.path.join(project_root, 'configs', 'config.yaml')
conf = OmegaConf.load(conf_path)

data_path =conf['outputs']['output_dir']

filename = f'adata_PC_with_ann_merged_v_{version}.h5ad' 
adata_cnv = sc.read_h5ad(os.path.join(data_path, filename))

adata_cnv.obs['pc_annotation'] = adata_cnv.obs['Populations'].copy().replace({'Normal_PC':'Healthy',
                                                                                'Interm':'Healthy_Like',
                                                                                'Normal_Pb':'Healthy_Like'})

adata_cnv.obs['Sample.Code.Cell'] = adata_cnv.obs['Sample.Code'].astype(str) + '.' + adata_cnv.obs['Populations'].astype(str)

# arch merge

path_save_8arch = '/home/projects/amit/annaku/repos/Blueprint/data/processed/nmf_outputs/renamed_with_8_separated/'
arch = pd.read_csv(os.path.join(path_save_8arch,f'arch_sample_v7_samplelevel_only_malignant_renamed_with_added_samples_v_{version}.csv'), index_col = 0)
arch[['Cluster', 'Cluster_exp']] = arch[['Cluster', 'Cluster_exp']].astype(str)

adata_cnv.obs['PID'] = 'z.' + adata_cnv.obs['Method'].astype(str) + '_' + adata_cnv.obs['Populations'].astype(str) + '_' + adata_cnv.obs["Sample.Code"].astype(str).str.lower()
adata_cnv.obs['PID'] = adata_cnv.obs['PID'].str.lower()
arch['PID'] = arch['PID'].str.lower() 
arch['patient'] = arch['patient'].str.lower()

adata_cnv.obs['index'] = adata_cnv.obs.index
adata_cnv.obs = pd.merge(adata_cnv.obs, arch[['Cluster_exp', 'Cluster', 'PID', 'prolif_high_0.1']], how='left', on='PID', validate='m:1')
adata_cnv.obs.index = adata_cnv.obs['index']

adata_cnv.obs = adata_cnv.obs.rename(columns = {'Cluster_exp':'arch', 'Cluster':'arch_with_8'})

arch_rename_dict = {'2.0':'MM1',
                    '3.0':'MM2',
                    '4.0':'MM3',
                    '6.0':'MM4',
                    '7.0':'MM5'}

adata_cnv.obs['arch'] = adata_cnv.obs['arch'].replace(arch_rename_dict)

adata_minibulk_combined = process_minibulk_cnv(adata_cnv,
                                            to_filter = False,
                                            reference_key = 'Populations',
                                        reference_cat=['Normal_PC'],
                                            preproc_standart = False,
                                            cat_to_group=['Sample.Code', 'Populations']
)

adata_minibulk_combined.obs = adata_minibulk_combined.obs.drop(columns = ['cnv_score', 'cnv_leiden'])

# separately for methods

savedir = data_path + '/infercnv_r_input/'

adata_minibulk_combined.obs["Populations_arch"] = (
    adata_minibulk_combined.obs["Populations"]
    + "_"
    + adata_minibulk_combined.obs["arch"].fillna("noarch")
)
adata_minibulk_combined.obs["Populations_arch_prolif"] = (
    adata_minibulk_combined.obs["Populations_arch"]
    + "_"
    + adata_minibulk_combined.obs["prolif_high_0.1"]
    .astype(str)
    .replace({"True": "Prolif", "False": "Non-Prolif"})
)

adata_minibulk_combined.obs['Populations_arch_prolif'].value_counts(dropna = False)

for method in ['MARS', 'SPID']:
    print(method)

    adata_minibulk_combined_method = adata_minibulk_combined[adata_minibulk_combined.obs['Method'] == method]

    adata_minibulk_combined_method_filtered = adata_minibulk_combined_method.copy()

    display(adata_minibulk_combined_method_filtered.obs['Populations_arch_prolif'].value_counts(dropna = False))

    ann = adata_minibulk_combined_method_filtered.obs[['Populations_arch_prolif', 'Sample.Code.Cell']]
    ann['Populations_arch_prolif'] = ann['Populations_arch_prolif'].replace({#'Healthy_noarch_nan':"Healthy",
                                                                #'Healthy_Like_noarch_nan':'Healthy_Like',
                                                                'Normal_PC_noarch_nan':'Normal_PC',
                                                                'Interm_noarch_nan':'Interm',
                                                                'Normal_Pb_noarch_nan':'Normal_Pb'
                                                                })
    display(ann['Populations_arch_prolif'].value_counts())

    ann.to_csv(savedir + f'ann_arch_prolif_{method}.txt', sep='\t', header=False)

    exp_df = pd.DataFrame(
        adata_minibulk_combined_method_filtered.X, 
        index=adata_minibulk_combined_method_filtered.obs.index,
        columns=adata_minibulk_combined_method_filtered.var.index
    ).T

    exp_df.to_csv(savedir + f'exp_arch_prolif_{method}.tsv', sep='\t')

    exp_df.head()

    adata_minibulk_combined_method_filtered.var = adata_minibulk_combined_method_filtered.var.merge(load_intervals(), left_index=True, right_index=True)
    display(adata_minibulk_combined_method_filtered.var.head())

    ann_gene = adata_minibulk_combined_method_filtered.var[['chromosome', 'start', 'end']]
    ann_gene.to_csv(savedir + 'ann_gene.txt', sep='\t', header=False)