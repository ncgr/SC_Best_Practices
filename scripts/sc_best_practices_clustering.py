import sys
import scanpy as sc


def run_clustering(h5ad_file, sample_id):
    """Clustering based on

       https://www.sc-best-practices.org/cellular_structure/clustering.html
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

    ## Neighbors for UMAP
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)

    ## Run leiden clustering
    sc.tl.leiden(adata)

    ## Leiden clustering at 3 resolutions controlling the KNN cluster density
    sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
    sc.pl.umap(
        adata,
        color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
        legend_loc="on data",
        save=f".{sample_id}.leiden_clustering.umap.png",
    )


if __name__ == "__main__":
    sample_id = sys.argv[1]  # get sample_id from command line, REMOVE LATER
    if not sample_id:
        sample_id = "2Cinoc_wholeroot"  # default for testing
    h5_input = f"./{sample_id}.dimensionality_reduction.h5ad"
    run_clustering(h5_input, sample_id)
