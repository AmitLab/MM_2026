import os
from typing import List
import sys
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import torch
from omegaconf import OmegaConf

project_root = Path(__file__).resolve().parent.parent 
sys.path.append(str(project_root / 'src'))

from io_utils import generate_path_in_output_dir
from train_scvi_model import generate_model_name, train_scvi_model, load_pp_adata_after_norm_and_hvg


def _get_adata_only_of_populations(adata: ad.AnnData, config, wanted_population_names: List[str]) -> ad.AnnData:
    assert config.annotation.cell_type_column in adata.obs.columns
    cells_mask = adata.obs[config.annotation.cell_type_column]. \
        apply(lambda x: x in wanted_population_names).replace(np.nan, False)
    assert sum(cells_mask) > 0, "no cells found of the given populations"
    return adata[cells_mask]


if __name__ == '__main__':
    torch.set_float32_matmul_precision("high")

    default_config = str(project_root / 'configs' / 'config.yaml')
    conf = OmegaConf.load(default_config)

    adata_for_training = load_pp_adata_after_norm_and_hvg(conf)

    population_names = ["Malignant", "Healthy"]
    adata_for_training = _get_adata_only_of_populations(adata_for_training, conf, population_names).copy()

    model = train_scvi_model(adata_for_training, batch_key=conf.scvi_settings.batch_key)

    model_name = generate_model_name(conf, extra_description=f"_only_{'_'.join(population_names)}")
    model_path = generate_path_in_output_dir(conf, model_name, add_date_timestamp=conf.outputs.add_timestamp)
    model.save(model_path, overwrite=True)
