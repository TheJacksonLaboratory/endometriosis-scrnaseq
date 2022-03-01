#!/usr/bin/env python
# -*- coding: utf-8 -*-
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def run_scvelo_paga(adata, n_pcs = 30, n_neighbors = 30, suffix="", subdir=""):
    adata_copy = adata.copy()
    
    print("Velocyto analysis for:")
    print(adata_copy.obs.subtypes.unique().to_list())
    
    #calculate moments
    scv.pp.neighbors(adata_copy, n_pcs=n_pcs, n_neighbors=n_pcs,  metric="correlation")
    scv.pp.moments(adata_copy, n_pcs=n_pcs, n_neighbors=n_neighbors, use_highly_variable=True)

    #Differenial kinetics analysis
    scv.tl.recover_dynamics(adata_copy)
    top_genes = adata_copy.var['fit_likelihood'].sort_values(ascending=False).index[:200]
    scv.tl.differential_kinetic_test(adata_copy, var_names=top_genes, groupby='subtypes')

    kwargs = dict(frameon=False, size=10, linewidth=1.5)
    scv.pl.scatter(adata_copy, basis=top_genes[:15], ncols=5, add_outline='fit_diff_kinetics', **kwargs)

    #velocity analysis
    scv.tl.velocity(adata_copy, mode="dynamical", diff_kinetics=True)
    scv.tl.velocity_graph(adata_copy,mode_neighbors="connectivities")
    scv.pl.velocity_embedding_stream(adata_copy, basis='umap', color=["subtypes"], legend_loc = "center right", title="")
    
    #get genes driving the stream in each cluster
    scv.tl.rank_velocity_genes(adata_copy, groupby="subtypes", min_corr=.3)
    
    ##PAGA
    scv.tl.paga(adata_copy, groups='subtypes')
    #save_adata(adata_copy, suffix=suffix, subdir=subdir)
    
    ## plot figure
    fig, ax = plt.subplots(1,1, figsize=(4,4))
    sc.pl.umap(adata_copy, color = ["subtypes"], alpha=.1,size=50, ax=ax, show=False)
    scv.pl.paga(adata_copy, basis='umap', size=50, alpha=.1, min_edge_width=1.5, node_size_scale=1, ax=ax)
    #fix_aspect_scatter_with_legend(fig)
    #save_figure(fig,pdir=pdir, name=name)
    
    #return fig
    #df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
    #df.style.background_gradient(cmap='Blues').format('{:.2g}')
    

