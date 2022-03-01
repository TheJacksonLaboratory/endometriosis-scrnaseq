#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Wedge, Circle
from matplotlib.legend_handler import HandlerPatch


class WedgeHandler(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        patch = Wedge(
            [xdescent, ydescent], width, 0, np.arctan(2*height/width)*180/np.pi, 
            facecolor=orig_handle[0], lw=0.9,
            edgecolor="k"
        )
        return [patch]


def _pie(ax, data, palette):
    # `normalize=True` is not needed in theory but 
    # good to handle floating point error
    ax.pie(
        data,
        normalize=True, 
        colors=palette,
        wedgeprops=dict(lw=1, ec="k"),
    )


def pies(p, palette=sns.mpl_palette("Dark2", 6)):
    celltypes = p["celltype"].unique()
    n_celltypes = len(celltypes)

    sample_types = p["sample_type"].unique()
    n_sample_types = len(sample_types)

    patients = p["patient"].unique()
    n_patients = len(patients)

    max_rows = (
        p.groupby(["sample_type", "patient"]).agg(sum) > 0
    ).groupby("sample_type").sum().max().values.ravel()[0]
    fig, axs = plt.subplots(
        max_rows + 1, n_sample_types, 
        dpi=200, figsize=(n_sample_types, n_patients+1),
        gridspec_kw=dict(hspace=0, wspace=0)
    )
    for ax in axs.flat:
        ax.set_axis_off()

    # Top row
    agg_data = p.groupby(["sample_type", "celltype"]).agg(np.mean)
    for k, st in enumerate(sample_types):
        idx = agg_data.index.get_level_values("sample_type") == st
        axs[0, k].set_axis_on()
        _pie(
            axs[0, k],
            agg_data.loc[idx].values.ravel(),
            palette
        )
        axs[0, k].set_title(st)

    # All other rows
    for j, st in enumerate(sample_types):
        i = 1
        for pat in patients:
            d = p[p.patient.isin([pat]) & p.sample_type.isin([st])]
            if sum(d["proportion"] > 0) < 1: 
                continue
            axs[i, j].set_axis_on()
            _pie(axs[i, j], d["proportion"], palette)
            axs[i, j].text(0.75, 0.85, pat, fontsize="x-small", transform=axs[i,j].transAxes)
            i += 1

    p0 = axs[0,0].get_position()
    p1 = axs[1,0].get_position()
    line_height = (p0.ymin + p1.ymax)/2

    line = plt.Line2D([0.05, 0.95],[line_height]*2, transform=fig.transFigure, color="black", lw=1)
    fig.add_artist(line)

    legend_elements = [(c, ct) for c, ct in zip(palette, celltypes)]
    axs[7, -2].legend(
        legend_elements, list(celltypes), ncol=1,
        loc="upper center", frameon=False, bbox_to_anchor=(1.25, 0.5),
        handler_map={tuple: WedgeHandler()}, handletextpad=0.5
    )

    return fig


def prep_props(adata, ct_key="celltype_main", pt_key="Patient_id", st_key="sample_type_rename"):
    n_ct = len(adata.obs[ct_key].unique())
    props = adata.obs.groupby([pt_key, st_key, ct_key]).size().reset_index()
    props.columns = ["patient", "sample_type", "celltype", "proportion"]
    norms = np.repeat(props.groupby(["patient", "sample_type"]).agg(sum).values.ravel(), n_ct)
    props.proportion /= norms
    return props
