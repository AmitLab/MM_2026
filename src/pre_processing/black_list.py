import logging
import re
from typing import List, Iterable

import anndata as ad


def get_mask_that_match_any_of_patterns(optional_strings: Iterable[str], patterns: List[str]) -> List[bool]:
    return [re.search(r"|".join(patterns), string) is not None for string in optional_strings]


def get_strings_that_match_no_pattern(optional_strings: Iterable[str], patterns: List[str]) -> List[str]:
    return [string for string in optional_strings if re.search(r"|".join(patterns), string) is None]


def drop_blacklist_genes_by_pattern(adata: ad.AnnData, blacklist_genes_patterns: List[str]) -> ad.AnnData:
    genes_to_keep = get_strings_that_match_no_pattern(adata.var_names, blacklist_genes_patterns)
    logging.info(f"droped {len(adata.var_names) - len(genes_to_keep)} genes "
                 f"that matched {blacklist_genes_patterns} patterns")
    return adata[:, genes_to_keep]
