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
