#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib

def fix_violin_edge_colors(ax, violin_edge_color="0.2", box_line_color="0.2", box_dot_color="white"):
    for c in ax.get_children():
        if isinstance(c, matplotlib.collections.PolyCollection):
            c.set_edgecolor("green")
        if isinstance(c, matplotlib.collections.PathCollection):
            c.set_color("blue")
        if isinstance(c, matplotlib.lines.Line2D):
            c.set_color("red")


def local_despine(ax, *, top=True, bottom=False, left=False, right=True):
    for side in ["top", "right", "left", "bottom"]:
        ax.spines[side].set_visible(not locals()[side])


def _remove_unused(s):
    s.cat = s.cat.remove_unused_categories()

def fix_main_adata_annos(adata):
    celltypes = [
        "epithelial", "stromal", "endothelial", "myeloid", "lymphocytes"
    ]
    sample_types = ["Ctrl", "EuE", "EcP", "EcPA", "EcO"]
    patients =  [f"C{k:02}" for k in range(1, 4) ] + [f"E{k:02}" for k in range(1, 12)]

    adata.obs["celltype_main"].cat.reorder_categories(celltypes, inplace=True)
    adata.obs["sample_type_rename"].cat.reorder_categories(sample_types, inplace=True)
    adata.obs["PID"].cat.reorder_categories(patients, inplace=True)

    _remove_unused(adata.obs["sample_type_rename"])
    _remove_unused(adata.obs["celltype_main"])
    _remove_unused(adata.obs["PID"])
