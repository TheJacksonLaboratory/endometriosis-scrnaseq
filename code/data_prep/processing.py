#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
import scanpy as sc

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import cdist
from scipy.stats import median_absolute_deviation
from scanpy.plotting._anndata import _prepare_dataframe


def merge_adata(
    *adatas,
    output_dir: str = os.getcwd(),
    sample_id: str = ""
):
    assert len(adatas) > 1
    first, *rest = adatas
    comb = first.concatenate(
        rest,
        join="outer",
        batch_key="library"
    )
    for key in ('10x_chemistry', 'sample_name', 'sampleid'): 
        comb.obs[key] = sum(
            ([adata.uns[key][0]] * len(adata) for adata in adatas), 
            []
        )
        
    comb.uns["sampleid"] = sample_id
    comb.uns["output_dir"] = output_dir
    tmp_var = pd.DataFrame(index=comb.var_names)
    for key in first.var_keys():
        dtype = first.var[key].dtype
        cols = comb.var.columns[comb.var.columns.str.startswith(key)]
        if dtype == bool:
            tmp_var[key] = comb.var[cols].fillna(False).any(axis=1)
        elif dtype == "object":
            tmp_var.loc[first.var_names, key] = first.var[key]
            for adata in rest:
                tmp_var.loc[adata.var_names, key] = adata.var[key]
    comb.var = tmp_var
    comb.var["total_counts"] = np.squeeze(np.asarray(comb.X.sum(axis=0)))
    comb.var["n_cells_by_counts"] = np.squeeze(np.asarray((comb.X > 0).sum(axis=0)))
    return comb


def preprocess(adata_raw, n_top_genes=1000, scale=False):
    #adata_raw.obs["n_counts_total"] = adata_raw.obs["n_counts"].copy()
    adata_raw.uns["raw_dtype"] = "normalized count"
    adata_raw.uns["n_top_genes"] = n_top_genes

    adata = adata_raw.copy()
    adata.layers["raw"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    adata.layers["normed"] = adata.X.copy()

    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="cell_ranger"
    )

    if scale:
        adata.layers["scaled"] = sc.pp.scale(adata.X)
        
    adata.var.loc[
        adata.var.cell_cycle | \
        adata.var.mitochondrial | \
        adata.var.ribosomal | \
        adata.var.hemoglobin | \
        adata.var.stress_response,
        "highly_variable"
    ] = False
    return adata


def load_pcs(adata, pcs):
    """
    This function serves to load in corrected PCs (e.g. from Harmony) and
    recompute the variance, variance ratio.  That is because the revised
    principal components may not maintain the same ordering (by variance
    explained) as the original PCs.  Hence, we want to load the new ones,
    recompute the variance explained, then (possibly reorder) the PCs so that
    the appropriate PCs are selected during downstream steps.

    This may have been addressed once Scanpy integrated Harmonypy, but that
    didn't exist at the time of starting this work.
    """
    new = adata.copy()
    new.obsm["X_pca_uncorrected"] = new.obsm["X_pca"].copy()
    total_var = np.var(new.X[:, new.var.highly_variable].toarray(), axis=0).sum()
    new.uns["pca"]["variance_unsorted"] = np.var(pcs, axis=0)
    new.uns["pca"]["variance_ratio_unsorted"] = np.var(pcs, axis=0) / total_var
    sort_order = np.argsort(new.uns["pca"]["variance_unsorted"])[::-1]
    new.obsm["X_pca"] = pcs[:, sort_order]
    new.uns["pca"]["variance"] = new.uns["pca"]["variance_unsorted"][sort_order]
    new.uns["pca"]["variance_ratio"] = new.uns["pca"]["variance_ratio_unsorted"][sort_order]
    return new

        
def detect_umap_doublets(adata, cluster_key, t=10):
    """
    No matter what, you sometimes get a cell or two that just doesn't localize
    with cells of the same cluster.  While these may not be true doublets, they
    are relatively few enough in number to be considered outliers and we want
    to remove them.  Rather than remove them via expression, quality metrics,
    etc. (as they remain after QC), just remove them if they are far enough
    away from the centroid of cluster density.
    """
    if "visible_doublet" in adata.obs_keys():
        adata.obs.drop("visible_doublet", inplace=True, axis=1)

    adata.obs["visible_doublet"] = False
    for cluster in adata.obs[cluster_key].cat.categories:
        inds = adata.obs_names[adata.obs[cluster_key] == cluster]
        coords = adata[inds].obsm["X_umap"]
        centroid = np.median(coords,axis=0)
        dists = cdist(coords, [centroid])
        thresh = t * median_absolute_deviation(dists)
        outliers = (dists > thresh).flatten()
        adata.obs.loc[inds[outliers], "visible_doublet"] = True


def reorder_clusters_hierarchical(adata, cluster_key, new_key="cluster"):
    """
    Reorder/rename clusters based on some sort of hierarchical similarity
    rather than the seemingly arbitrary labels supplied by graph-based methods.
    Also start clustering at 1 rather than 0.

    Note this only changes labels, and does not "recluster" cells.
    """
    cats, tidy = _prepare_dataframe(
        adata,
        adata.var_names[adata.var.highly_variable],
        groupby=cluster_key,
        use_raw=False
    )
    mean_obs = tidy.groupby(level=0).mean()

    Z = linkage(mean_obs, method="ward", metric="euclidean", optimal_ordering=True)
    cluster_order = leaves_list(Z)
    new_clusters = (np.arange(len(cats)) + 1).astype(str)

    adata.obs[new_key] = adata.obs[cluster_key].map(
        dict(zip(cats[cluster_order], new_clusters))
    ).astype("category")
    adata.obs[new_key].cat.reorder_categories(new_clusters, inplace=True)

    return adata
