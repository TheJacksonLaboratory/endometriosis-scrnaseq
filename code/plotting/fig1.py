#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from plotting.chord import sankey


def set_thousands_ticks(ax, ticks, axis="x"):
    _ax = ax.xaxis
    f = ax.set_xlim
    if axis != "x":
        _ax = ax.yaxis
        f = ax.set_ylim
    _ax.set_ticks(ticks)
    labs = []
    for k in ticks:
        l = k
        m = k%1000
        if k >= 1000:
            if m > 0:
                l = f"{k/1000:.1f}k"
            else:
                l = f"{k//1000}k"
        labs.append(l)
    _ax.set_ticklabels(labs)
    f(ticks[0]*0.99, ticks[-1]*1.01)


def summary_figure_1d(
    adata, 
    left_margin, 
    right_margin,
    *,
    x_palette=[],
    left_palette=[],
    right_palette=[],
    cmap="magma",
    use_log=False,
    violin_limits={}
):
    plot_limits = {
        "patient": {"umis": (1000, 10000, 50000), "genes": (500, 3000, 10000), "cells": (0, 5000, 10000)},
        "sample_type": {"umis": (1000, 10000, 50000), "genes": (500, 3000, 10000), "cells": (0, 10000, 20000, 30000)},
    }
    plot_limits.update(**violin_limits)


    width = 3#len(adata.obs[x_column].unique())
    height = len(adata.obs[left_margin].unique())
    fig = plt.figure(facecolor="white", dpi=300, figsize=(8, 6))
    gs = plt.GridSpec(
        5, 9, 
        width_ratios=[2,2,2,0.5,6,0.5,2,2,2], 
        height_ratios=[3,2*height,2,2,1], 
        hspace=0.1, wspace=0.2
    )
    
    plt.rcParams["axes.axisbelow"] = True
    
    pat_counts = fig.add_subplot(gs[1,2])
    pat_labels = fig.add_subplot(gs[1,3], sharey=pat_counts)
    pat_umis   = fig.add_subplot(gs[1,0], sharey=pat_counts)
    pat_genes  = fig.add_subplot(gs[1,1], sharey=pat_counts)

    explainer_ax = fig.add_subplot(gs[0, 4])
    
    sankey_ax = fig.add_subplot(gs[1,4])
    st_umis   = fig.add_subplot(gs[1,8])
    st_genes  = fig.add_subplot(gs[1,7], sharey=st_umis)
    st_labels = fig.add_subplot(gs[1,5], sharey=st_umis)
    st_counts = fig.add_subplot(gs[1,6], sharey=st_umis)

    genes_axs = [pat_genes, st_genes]
    umis_axs = [pat_umis, st_umis]
    count_axs = [pat_counts, st_counts]
    for a in genes_axs + umis_axs + count_axs:
        a.grid(zorder=-100, lw=0.25, color="0.5")

    for col in [left_margin, right_margin]:
        adata.obs[col] = adata.obs[col].astype("category")
    
    pats = adata.obs[left_margin].cat.categories
    sts = adata.obs[right_margin].cat.categories
    
    count_kws = dict(lw=0.5, edgecolor="0.2")
    violin_kws = dict(inner="box", linewidth=0.7, linecolor="0.2", scale="width", cut=0)
    
    ## PATIENT SPECIFIC STUFF
    sns.countplot(data=adata.obs, y=left_margin, ax=pat_counts, palette=left_palette, **count_kws)
    
    for k, pat in enumerate(pats):
        pat_labels.text(0.5, k, pat, va="center", ha="center")
    
    sns.violinplot(data=adata.obs, y=left_margin, x="total_counts", ax=pat_umis, palette=left_palette, **violin_kws)
    sns.violinplot(data=adata.obs, y=left_margin, x="n_genes_by_counts", ax=pat_genes, palette=left_palette, **violin_kws)

    ## SANKY
    patients_to_sample_types = pd.pivot_table(
        adata.obs, index=right_margin, columns=left_margin, fill_value=0, aggfunc="size"
    ).astype(int).T
    sankey(
        sankey_ax, 
        patients_to_sample_types, 
        split_dist=0.90, rect_thickness=0.1,
        angle1=-90, angle2=-90,
        row_palette=left_palette,
        col_palette=right_palette,
        add_labels=False,
        add_edges=True
    )
    sankey_ax.set_xlim(-1, 1)
    sankey_ax.set_ylim(-1, 0)
    sankey_ax.set_clip_on(False)
    
    ## SAMPLE TYPE SPECIFIC STUFF  
    sns.countplot(data=adata.obs, y=right_margin, ax=st_counts, palette=right_palette, **count_kws)
    st_counts.invert_yaxis() # bottom bar is the first entry
    
    for k, st in enumerate(sts):
        st_labels.text(0.5, k, st, va="center", ha="center")
        
    sns.violinplot(data=adata.obs, y=right_margin, x="total_counts", ax=st_umis, palette=right_palette, **violin_kws)
    sns.violinplot(data=adata.obs, y=right_margin, x="n_genes_by_counts", ax=st_genes, palette=right_palette, **violin_kws)
    
    ## ALL LABEL, TICK, ETC ADJUSTMENTS
    for ax in fig.get_axes():
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    st_counts.yaxis.set_label_position("right")
    st_genes.yaxis.set_label_position("right")
    st_umis.yaxis.set_label_position("right")

    if use_log:
        st_umis.set_xscale("log")
        st_genes.set_xscale("log")
        pat_umis.set_xscale("log")
        pat_genes.set_xscale("log")
    
    set_thousands_ticks(st_counts, plot_limits["sample_type"]["cells"])
    set_thousands_ticks(st_genes, plot_limits["sample_type"]["genes"])
    set_thousands_ticks(st_umis, plot_limits["sample_type"]["umis"])
    set_thousands_ticks(pat_counts, plot_limits["patient"]["cells"])
    set_thousands_ticks(pat_genes, plot_limits["patient"]["genes"])
    set_thousands_ticks(pat_umis, plot_limits["patient"]["umis"])
    
    pat_umis.invert_xaxis()
    pat_genes.invert_xaxis()
    pat_counts.invert_xaxis()

    pat_labels.tick_params(bottom=False, left=False, right=False, top=False)
    st_labels.tick_params(bottom=False, left=False, right=False, top=False)
    explainer_ax.tick_params(bottom=False, left=False, right=False, top=False)

    pat_counts.tick_params(bottom=True, left=False, right=False, top=False)
    st_counts.tick_params(bottom=True, left=False, right=False, top=False)
    pat_genes.tick_params(bottom=True, left=False, right=False, top=False)
    st_genes.tick_params(bottom=True, left=False, right=False, top=False)
    pat_umis.tick_params(bottom=True, left=False, right=False, top=False)
    st_umis.tick_params(bottom=True, left=False, right=False, top=False)
    
    for ax in fig.get_axes():
        ax.tick_params(pad=1, labelsize=6)

    sns.despine(ax=pat_labels, top=True, bottom=True, left=True, right=True)
    sns.despine(ax=pat_counts, left=True, bottom=False, top=True, right=False, trim=True)
    sns.despine(ax=pat_umis, top=True, left=True, right=True, trim=True)
    sns.despine(ax=pat_genes, top=True, left=True, right=True, trim=True)

    sns.despine(ax=st_labels, top=True, bottom=True, right=True, left=True)
    sns.despine(ax=st_counts, right=True, top=True, bottom=False, left=False, trim=True)
    sns.despine(ax=st_umis, top=True, left=False, right=True, trim=True)
    sns.despine(ax=st_genes, top=True, left=True, right=True, trim=True)
    
    sns.despine(ax=explainer_ax, top=True, left=True, right=True, trim=True)
    explainer_ax.text(0.5, 0.25, r"Patient $\rightarrow$ Sample Type", ha="center", transform=explainer_ax.transAxes)
    
    pat_umis.set_ylabel("Patient", size=10)
    #pat_labels.set_title("PID", size=10)
    st_umis.set_ylabel("Sample Type", size=10, rotation=270, labelpad=10)

    pat_counts.set_xlabel("Cells", size=8)
    pat_genes.set_xlabel("Genes", size=8)
    pat_umis.set_xlabel("UMIs", size=8)
    st_genes.set_xlabel("Genes", size=8)
    st_umis.set_xlabel("UMIs", size=8)
    st_counts.set_xlabel("Cells", size=8)

    return fig


def summary_figure_1e(
    adata, 
    x_column, 
    left_margin, 
    right_margin,
    *,
    x_palette=[],
    left_palette=[],
    right_palette=[],
    cmap="magma",
    use_log=False,
):
    width = len(adata.obs[x_column].unique())
    height = len(adata.obs[left_margin].unique())
    fig = plt.figure(facecolor="white", dpi=300, figsize=(8, 4))
    gs = plt.GridSpec(
        3, 6, 
        width_ratios=[3,2,width,1.5,1.5,1], 
        height_ratios=[height, height, 1], 
        hspace=0, wspace=0
    )
    
    plt.rcParams["axes.axisbelow"] = True
    
    ct_main = fig.add_subplot(gs[:2,2])
    ct_counts = fig.add_subplot(gs[:2,0], sharey=ct_main)
    ct_labels = fig.add_subplot(gs[:2,1], sharey=ct_main)
    
    st_labels = fig.add_subplot(gs[2,2], sharex=ct_main)
    
    sankey_ax = fig.add_subplot(gs[:2,3:5])
    imc_labels = fig.add_subplot(gs[0,5])
    ct_cbar = fig.add_subplot(gs[2, 4:])

    for a in [ct_counts,]:
        a.grid(zorder=-100, lw=0.5, color="0.5")

    for col in [x_column, left_margin]:
        adata.obs[col] = adata.obs[col].astype("category")
    
    sts = adata.obs[x_column].cat.categories
    cts = adata.obs[left_margin].cat.categories
    imcs = right_margin
    
    ## MAIN PLOT
    xs = np.arange(width)
    main_data = pd.pivot_table(adata.obs, columns=x_column, index=left_margin, aggfunc="size")
    main_pat_norm = (main_data / main_data.sum(axis=0) * 100).astype(int)
    main_pat_max = main_pat_norm.values.max()
    
    im = ct_main.imshow(
        main_pat_norm.values,
        #aspect="auto", 
        origin="upper",
        cmap=cmap
    )
    fig.colorbar(im, cax=ct_cbar, orientation="horizontal")
    
    for i, r in enumerate(main_pat_norm.index):
        for j, c in enumerate(main_pat_norm.columns):
            color = im.cmap(main_pat_norm.loc[r,c] / main_pat_max)
            lum = sns.utils.relative_luminance(color)
            text_color = ".15" if lum > .408 else "w"
            ct_main.text(j, i, (main_pat_norm.loc[r,c]), color=text_color, fontsize=8, va="center", ha="center")


    ## CELLTYPE SPECIFIC STUFF
    count_kws = dict(lw=0.5, edgecolor="0.2")
    sns.countplot(data=adata.obs, y=left_margin, ax=ct_counts, palette=left_palette, **count_kws)
    
    for k, ct in enumerate(cts):
        ct_labels.text(0.5, k, ct, va="center", ha="center")

    ## ST SPECIFIC STUFF
    for k, st in enumerate(sts):
        st_labels.text(k, 0.5, st, va="center", ha="center", rotation=0, rotation_mode="anchor")
    
    # IMC STUFF
    if len(right_palette) < len(right_margin):
        right_palette = sns.mpl_palette("Set1", len(right_margin))
    for k, (key, color) in enumerate(zip(right_margin, right_palette)):
        imc_labels.text(0.1, k, key, color=color, va="center", ha="left", )
    imc_labels.set_ylim(-0.5, len(right_margin)-0.5)
    imc_labels.invert_yaxis()
    
    ## SANKY
    imc_panel = pd.DataFrame(
        np.eye(len(right_margin)),
        columns = right_margin,
        index=cts
    )
    sankey(
        sankey_ax, 
        imc_panel, 
        split_dist=0.90, rect_thickness=0.1,
        angle1=-90, angle2=-90,
        row_palette=left_palette,
        col_palette=right_palette,
        add_labels=False,
        add_edges=True,
        norm = True,
        receiver_length=0.5
    )
    sankey_ax.set_xlim(-1, 1)
    sankey_ax.set_ylim(-1, 0)
    sankey_ax.set_clip_on(False)
    
    ## ALL LABEL, TICK, ETC ADJUSTMENTS
    for ax in fig.get_axes():
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    ct_counts.invert_xaxis()
    #ct_counts.invert_yaxis()
    ct_counts.set_xticks([0, 15000, 30000, 45000])
    ct_counts.set_xticklabels(["0", "15k", "30k", "45k"])
    ct_counts.yaxis.set_label_position("right")
    ct_counts.yaxis.tick_right()
    
    for ax in fig.get_axes():
        ax.tick_params(pad=1, labelsize=6)
    
    all_true = dict(top=True, right=True, bottom=True, left=True)
    all_false = dict(top=False, right=False, bottom=False, left=False)

    ct_main.tick_params(bottom=False, left=False)
    ct_labels.tick_params(**all_false)
    st_labels.tick_params(**all_false)
    imc_labels.tick_params(**all_false)
    sankey_ax.tick_params(**all_false)

    sns.despine(ax=ct_main, bottom=True, left=True)
    sns.despine(ax=ct_counts, top=True, left=True, right=False, trim=True)
    sns.despine(ax=ct_labels, **all_true)
    sns.despine(ax=st_labels, **all_true)
    sns.despine(ax=imc_labels, **all_true)
    sns.despine(ax=sankey_ax, **all_true)
    
    ct_counts.set_xlabel("Cells", size=8)
    ct_cbar.set_xlabel("Relative proportion within each sample type", size=8)

    return fig
