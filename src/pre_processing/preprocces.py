import re

import anndata as ad


def drop_bad_genes(adata: ad.AnnData) -> ad.AnnData:
    bad_genes_regexs = [r"^Rps", r"^Rpl", r".*Rik$", r"^Gm[0-9]", r"^Gm[0-9]", r"^AC[0-9][0-9]", r"^mt-", r"^Mtmr",
                        r"^Mtnd", "Neat1", "Malat1"]
    adata = adata[:, [gn for gn in adata.var_names if re.match("|".join(bad_genes_regexs), gn)]]
    return adata
