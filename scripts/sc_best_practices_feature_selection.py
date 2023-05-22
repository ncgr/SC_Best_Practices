import sys
import scanpy as sc
import anndata2ri
import logging
import numpy as np

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro


def deviance_feature_selection(adata):
    """Calls feature selection with deviance from non-normalized counts"""
    ro.r("""
        library(scry)
        sce = devianceFeatureSelection(adata, assay="X")
        """)
    binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T

    ## get top 4000 highly deviant genes
    idx = binomial_deviance.argsort()[-4000:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    ## set results as variables in adata
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance


def run_feature_selection(h5ad_file, sample_id):
    """Feature selection based on 
    
       https://www.sc-best-practices.org/preprocessing_visualization/feature_selection.html
    """
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(
        dpi=80,
        facecolor="white",
        frameon=False,
    )

    rcb.logger.setLevel(logging.ERROR)
    ro.pandas2ri.activate()
    anndata2ri.activate()

    adata = sc.read(
        filename=h5ad_file,
    )
    ro.globalenv["adata"] = adata

    ## Deviance feature selection
    deviance_feature_selection(adata)

    ## Write output h5ad file
    adata.write(f"./{sample_id}.feature_selection.h5ad")


if __name__ == "__main__":
    sample_id = sys.argv[1]  # get sample_id from command line, REMOVE LATER
    if not sample_id:
        sample_id = "2Cinoc_wholeroot"  # default for testing
    h5_input = f"./{sample_id}.normalization.h5ad"
    run_feature_selection(h5_input, sample_id)
