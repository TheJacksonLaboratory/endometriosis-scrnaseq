#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scanpy.plotting._anndata import _prepare_dataframe

def report_in_out_group_fractions(
    adata: sc.AnnData,
    key=None,
    groupby=None,
    use_raw=True,
    log=True,
) -> None:
    if key is None:
        key = 'rank_genes_groups'
    if groupby is None:
        groupby = str(adata.uns[key]['params']['groupby'])
    gene_names = pd.DataFrame(adata.uns[key]['names'])

    fraction_in_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )

    fraction_out_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )

    for cluster in gene_names.columns:
        var_names = gene_names[cluster].values
        adata.obs['__is_in_cluster__'] = pd.Categorical(adata.obs[groupby] == cluster)
        categories, obs_tidy = _prepare_dataframe(
            adata,
            var_names,
            groupby='__is_in_cluster__',
            use_raw=use_raw,
        )
        mean_obs = obs_tidy.groupby(level=0).mean()
        obs_bool = obs_tidy.astype(bool)
        fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()
        fraction_in_cluster_matrix.loc[:, cluster] = fraction_obs.loc[True].values
        fraction_out_cluster_matrix.loc[:, cluster] = fraction_obs.loc[False].values
    adata.obs.drop(columns='__is_in_cluster__')

    return fraction_in_cluster_matrix, fraction_out_cluster_matrix

_helper_outer = lambda df, name: pd.melt(
    df, value_name=name, var_name="cluster"
).set_index("cluster")
_helper = lambda adata, col, name: _helper_outer(
    pd.DataFrame.from_records(adata.uns["rank_genes_groups"][col]), name
)

def compute_marker_genes(
    adata, outdir, groupby="cluster", use_raw=True, method="wilcoxon",
    suffix="", save=True
):
    sc.tl.rank_genes_groups(adata, groupby=groupby, use_raw=use_raw, method="wilcoxon", n_genes=250)
    in_groups, out_groups = report_in_out_group_fractions(adata, groupby=groupby)
    sizes = adata.obs.groupby(groupby).size()
    markers = pd.concat(
        (
            _helper(adata, "names", "gene"),
            _helper(adata, "scores", "score"),
            _helper(adata, "logfoldchanges", "logfc"),
            _helper(adata, "pvals", "pval"),
            _helper(adata, "pvals_adj", "pval_fdr"),
            _helper_outer(in_groups, "in_group_fraction"),
            _helper_outer(out_groups, "out_group_fraction"),
            _helper_outer(out_groups < 0.2, "specific_20%"),
            _helper_outer(out_groups < 0.05, "specific_5%"),
            _helper_outer(out_groups < 0.01, "specific_1%")
        ),
        axis=1
    )
    markers["n_cells"] = sizes.loc[markers.index]
    # remove ribosomal genes
    markers = markers.loc[~markers.gene.str.match("^M?RP[LS]"), :]
    
    if save:
        suffix = f"-{suffix}" if suffix else ""
        #outpath = Path(ANALYSIS_DIR / "markers" / f"{sample_id}{suffix}-nn.csv")
        outpath = f"{outdir}/markers/{adata.uns['sampleid']}{suffix}-markers.csv"
        markers.to_csv(outpath, index=True, header=True)
    else:
        return markers
