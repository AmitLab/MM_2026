from pathlib import Path

from omegaconf import DictConfig, OmegaConf

from data_loading.anndata_loaders import FromPlatesDataLoader, AnnDataLoader, MultiMethodFromPlatesDataLoader


# def create_dataloader_from_config(config: DictConfig) -> AnnDataLoader:
#     if OmegaConf.is_list(config.data_loading.sc_sequencing.sc_sequencing_data_dir):
#         sc_dirs = [Path(path) for path in config.data_loading.sc_sequencing.sc_sequencing_data_dir]
#         plate_data_paths = [Path(path) for path in config.data_loading.plates.plates_data_path]
#         return MultiMethodFromPlatesDataLoader(
#             sc_data_dirs=sc_dirs,
#             plates_data_paths=plate_data_paths,
#             plate_id_col_name=config.data_loading.plates.plate_id_column_name,
#             plates_data_transform_functions=[])
#     return FromPlatesDataLoader(sc_data_dir=Path(config.data_loading.sc_sequencing.sc_sequencing_data_dir),
#                                 plates_data_path=Path(config.data_loading.plates.plates_data_path),
#                                 plate_id_col_name=config.data_loading.plates.plate_id_column_name,
#                                 plates_data_transform_functions=[])

def create_dataloader_from_config(config: DictConfig) -> AnnDataLoader:
    sc_data_dir = config.data_loading.sc_sequencing.sc_sequencing_data_dir
    plates_data_path = config.data_loading.plates.plates_data_path
    # universal check for list-like objects
    is_list_like = hasattr(sc_data_dir, '__iter__') and not isinstance(sc_data_dir, str)
    
    if is_list_like:
        sc_dirs = [Path(path) for path in sc_data_dir]
        plate_data_paths = [Path(path) for path in config.data_loading.plates.plates_data_path]
        return MultiMethodFromPlatesDataLoader(
            sc_data_dirs=sc_dirs,
            plates_data_paths=plate_data_paths,
            plate_id_col_name=config.data_loading.plates.plate_id_column_name,
            plates_data_transform_functions=[])
    return FromPlatesDataLoader(sc_data_dir=Path(config.data_loading.sc_sequencing.sc_sequencing_data_dir),
                                plates_data_path=Path(config.data_loading.plates.plates_data_path),
                                plate_id_col_name=config.data_loading.plates.plate_id_column_name,
                                plates_data_transform_functions=[])
