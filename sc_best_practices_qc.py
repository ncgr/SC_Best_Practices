import sys
import logging
import anndata2ri
import numpy as np
import scanpy as sc
import seaborn as sns
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from scipy.stats import median_abs_deviation


def run_doublet_detection():
    """Uses scDblFinder to find the number of doublets"""
    ro.r("""
        library(Seurat)
        library(scater)
        library(scDblFinder)
        library(BiocParallel)

        set.seed(123)
        sce = scDblFinder(
            SingleCellExperiment(
                list(counts=data_mat),
            )
        )
        doublet_score = sce$scDblFinder.score
        doublet_class = sce$scDblFinder.class
     """)


def run_soupx():
    """Uses rpy2.robjects to run soupX"""
    ro.r("""
        library(SoupX)
        set.seed(123)
        # specify row and column names of data
        rownames(data) = genes
        colnames(data) = cells
        # ensure correct sparse format for table of counts and table of droplets
        data <- as(data, "sparseMatrix")
        data_tod <- as(data_tod, "sparseMatrix")
        
        # Generate SoupChannel Object for SoupX
        sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)
        
        # Add extra meta data to the SoupChannel object
        soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
        sc = setSoupProfile(sc, soupProf)
        # Set cluster information in SoupChannel
        sc = setClusters(sc, soupx_groups)
        
        # Estimate contamination fraction
        sc  = autoEstCont(sc, doPlot=FALSE)
        # Infer corrected table of counts and rount to integer
        out = adjustCounts(sc, roundToInt = TRUE)
    """)


def rna_correction(adata, adata_raw):
    """Ambient RNA correction method"""
    rcb.logger.setLevel(logging.ERROR)
    ro.pandas2ri.activate()
    anndata2ri.activate()
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")

    ## Preprocess variables for SoupX
    soupx_groups = adata_pp.obs["soupx_groups"]
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    ## Get Rawdata for SoupX
    adata_raw.var_names_make_unique()
    data_tod = adata_raw.X.T

    ## set robjects for R inputs
    ro.globalenv["genes"] = genes
    ro.globalenv["cells"] = cells
    ro.globalenv["data"] = data
    ro.globalenv["data_tod"] = data_tod
    ro.globalenv["soupx_groups"] = soupx_groups


def is_outlier(adata, metric, nmads):
    """Decides if outlier based on MAD"""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def run_qc(h5_file, h5_file_raw, sample_id):
    """Takes a matrix.h5 file and processes it through QC based on

       https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
    """
    ## setup
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(
        dpi=80,
        facecolor="white",
        frameon=False,
    )
    ## read file
    adata = sc.read_10x_h5(
        filename=h5_file,
    )
    adata_raw = sc.read_10x_h5(
        filename=h5_file_raw,
    )
#    print(adata)
    ## make names unique
    adata.var_names_make_unique()
    adata_raw.var_names_make_unique()
#    print(adata)
    
    ## mitochondrial for A17
    adata.var["mt"] = adata.var_names.str.contains("A17MT")
    ## chloroplast for A17
    adata.var["cp"] = adata.var_names.str.contains("A17CP")
#    print(adata)
    
    ## scanpy QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "cp"], inplace=True, percent_top=[20], log1p=True
    )
    print(adata)  # print data in final form
    
    ## Plots
    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    # sc.pl.violin(adata, 'total_counts')
    p2 = sc.pl.violin(adata, "pct_counts_mt", save=f".{sample_id}.pct_counts_mt.png")
    p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save=f".{sample_id}.total_countsXn_genes_by_count.png")
    p1.savefig(f"./figures/{sample_id}.total_counts.png")
    ## outliers 3MADs MAD5 threshold
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    ## Cells with MT >= 8% Removed
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    print(adata.obs.outlier.value_counts())
    print(f"\nTotal number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save=f".{sample_id}.total_countsXn_genes_by_count.final.png")

    ## Run Ambient RNA correction
    rna_correction(adata, adata_raw)

    ## Run soupX in R using rpy2
    run_soupx()

    ## Add SoupX output to adata
    adata.layers["counts"] = adata.X
    adata.layers["soupX_counts"] = ro.globalenv["out"].T
    adata.X = adata.layers["soupX_counts"]

    print(f"\nTotal number of genes: {adata.n_vars}")

    ## Min 20 cells - filters out 0 count genes
    sc.pp.filter_genes(adata, min_cells=20)
    print(f"Number of genes after cell filter: {adata.n_vars}")

    ## Doublet Detection
    data_mat = adata.X.T
    ro.globalenv["data_mat"] = data_mat
    run_doublet_detection()
    adata.obs["scDblFinder_score"] = ro.globalenv["doublet_score"]
    adata.obs["scDblFinder_class"] = ro.globalenv["doublet_class"]
    print(adata.obs.scDblFinder_class.value_counts())
    
    ## output
    adata.write(f"./{sample_id}.qc.h5ad")


if __name__ == "__main__":
    sample_id = sys.argv[1]  # get sample_id from command line, REMOVE LATER
    if not sample_id:
        sample_id = "2Cinoc_wholeroot"  # default for testing
    h5_input = f"./{sample_id}/outs/filtered_feature_bc_matrix.h5"
    h5_input_raw = f"./{sample_id}/outs/raw_feature_bc_matrix.h5"
    run_qc(h5_input, h5_input_raw, sample_id)
