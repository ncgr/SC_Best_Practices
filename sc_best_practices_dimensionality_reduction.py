import sys
import scanpy as sc


def run_dimensionality_reduction(h5ad_file, sample_id, normalization="log1p_norm"):
    """PCA, t-SNE, UMAP. 

       https://www.sc-best-practices.org/preprocessing_visualization/dimensionality_reduction.html
    """
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(
        dpi=80,
        facecolor="white",
        frameon=False,
    )

    adata = sc.read(
        filename=h5ad_file,
    )

    ## for now use log1p_norm as defailt we also have scran and pearson residuals
    adata.X = adata.layers[normalization]

    ## setting highly variable as highly deviant to use scanpy 'use_highly_variable' argument in sc.pp.pca
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)  # PCA
    sc.pl.pca_scatter(adata, color="total_counts", save=f".{sample_id}.{normalization}.pca_scatter.png")

    ## t-SNE
    sc.tl.tsne(adata, use_rep="X_pca")
    sc.pl.tsne(adata, color="total_counts", save=f".{sample_id}.{normalization}.tsne.png")

    ## UMAP uses Neighbor graph
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="total_counts", save=f".{sample_id}.{normalization}.umap.png")


if __name__ == "__main__":
    sample_id = sys.argv[1]  # get sample_id from command line, REMOVE LATER
    if not sample_id:
        sample_id = "2Cinoc_wholeroot"  # default for testing
    h5_input = f"./{sample_id}.feature_selection.h5ad"
    run_dimensionality_reduction(h5_input, sample_id, "log1p_norm")
