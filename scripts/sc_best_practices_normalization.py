import sys
import logging
import scanpy as sc
import numpy as np
import seaborn as sns
import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from matplotlib import pyplot as plt
from scipy.sparse import issparse, csr_matrix



def shifted_logarithm(adata, sample_id):
    """Applies the fast normalization shifted logarithm
    
       This is based on the delta method.

       Uses the scanpy default L = median raw count depth
    """
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)

    ## log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
    axes[1].set_title("Shifted logarithm")
    plt.savefig(f"./figures/{sample_id}.shifted_logarithm_normalization.png")


def scran_size_factor_estimation(adata, sample_id):
    """Scran’s pooling-based size factor estimation method
    
       This is based on the delta method.
    """

    ## Preliminary clustering for differentiated normalisation
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups")

    data_mat = adata_pp.X.T
    ## convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
    if issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
    ro.globalenv["data_mat"] = data_mat
    ro.globalenv["input_groups"] = adata_pp.obs["groups"]

    ## Size factros for cell groups
    ro.r("""
            library(scran)
            library(BiocParallel)
            size_factors = sizeFactors(
                computeSumFactors(
                    SingleCellExperiment(
                        list(counts=data_mat)), 
                        clusters = input_groups,
                        min.mean = 0.1,
                        BPPARAM = MulticoreParam()
                )
            )
        """)
    adata.obs["size_factors"] = ro.globalenv["size_factors"]
    scran = adata.X / adata.obs["size_factors"].values[:, None]
    adata.layers["scran_normalization"] = csr_matrix(sc.pp.log1p(scran))  # set layer for scran normalization
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(
        adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
    )
    axes[1].set_title("log1p with Scran estimated size factors")
    plt.savefig(f"./figures/{sample_id}.log1p_scran_normalization.png")


def pearson_residuals(adata, sample_id):
    """https://www.sc-best-practices.org/preprocessing_visualization/normalization.html#analytic-pearson-residuals"""
    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
    adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(
        adata.layers["analytic_pearson_residuals"].sum(1), bins=100, kde=False, ax=axes[1]
    )
    axes[1].set_title("Analytic Pearson residuals")
    plt.savefig(f"./figures/{sample_id}.pearson_residuals_normalization.png")


def run_normalization(h5ad_file, sample_id):
    """Normalize QC h5ad output file"""
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

    ## Plot histogram of total_counts from adata
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False)
    p1.figure.savefig(f"./figures/{sample_id}.total_counts.normalization_input.png")

    ## shifted logarithm
    print("Shifted Logarithm Normalization")
    shifted_logarithm(adata, sample_id)  # delta method

    ## scran size factor estimation
    print("Log1p Scran Normalization")
    scran_size_factor_estimation(adata, sample_id)  # delta method

    ## Pearson residuals
    print("Pearson Residuals Normalization")
    pearson_residuals(adata, sample_id)  # scanpy
    adata.write(f"{sample_id}.normalization.h5ad")
    for k in adata.layers:
        print(k)


if __name__ == "__main__":
    sample_id = sys.argv[1]  # get sample_id from command line, REMOVE LATER
    if not sample_id:
        sample_id = "2Cinoc_wholeroot"  # default for testing
    h5_input = f"./{sample_id}.qc.h5ad"
    run_normalization(h5_input, sample_id)
