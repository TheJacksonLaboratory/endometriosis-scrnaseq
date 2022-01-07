#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.path import Path as mpl_path
import matplotlib.patches as patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D

import numpy as np
import pandas as pd


def build_flux_mat(df, senders=None, receivers=None, piv_kws={"aggfunc":"size"}, dtype=int):
    flux = df.pivot_table(index="celltype_a", columns="celltype_b", fill_value=0, **piv_kws).astype(int)

    possible_cts = flux.columns.union(flux.index)
    missing_cols = possible_cts.difference(flux.columns)
    missing_rows = possible_cts.difference(flux.index)
    for col in missing_cols:
        flux[col] = 0
    for row in missing_rows:
        flux.loc[row, :] = 0

    flux = flux.reindex(flux.columns)

    if senders is not None:
        flux.loc[~flux.index.isin(senders), :] = 0
    if receivers is not None:
        flux.loc[:, ~flux.columns.isin(receivers)] = 0

    return flux.astype(dtype)


def rot_mat(a):
    return np.array([np.cos(a), -np.sin(a), np.sin(a), np.cos(a)]).reshape(2, 2)


def rectangular_ideogram(
    ax,
    start,
    end,
    thickness=0.2,
    angle=0,
    origin_offset=0,
    facecolor=None,
    edgecolor=None,
    label="",
    label_ha="left"
):
    start, end = start, end
    xy = (start, 0)
    width = end - start
    height = thickness
    # Rectange(angle=angle) is about xy, so we need to compute offsets ourself
    transform = Affine2D().rotate_deg_around(0, 0, angle)             + Affine2D().translate(origin_offset, 0)             + ax.transData
    patch = patches.Rectangle(xy, width, height, facecolor=facecolor, edgecolor=edgecolor,
                             transform=transform)
    ax.add_patch(patch)
    ax.text((end+start)/2, thickness, s=label, transform=transform, fontsize=6, ha=label_ha)


def chord(
    ax,
    start1=-1, 
    end1=-2, 
    start2=1, 
    end2=2,
    angle1=0,
    angle2=0,
    x_offset=0,
    chordwidth=0.7, 
    facecolor=None,
    edgecolor=None
):
    assert start1 <= 0
    assert end1 <= 0
    
    angle1 *= np.pi / 180
    angle2 *= np.pi / 180
    x1, y1 = abs(start1) * np.cos(angle1) - x_offset, abs(start1) * np.sin(angle1)
    x2, y2 = abs(start2) * np.cos(angle2) + x_offset, abs(start2) * np.sin(angle2)
    x3, y3 = abs(end2) * np.cos(angle2) + x_offset, abs(end2) * np.sin(angle2)
    x4, y4 = abs(end1) * np.cos(angle1) - x_offset, abs(end1) * np.sin(angle1)

    codes = [
        mpl_path.MOVETO,
        mpl_path.CURVE4,
        mpl_path.CURVE4,
        mpl_path.CURVE4,
        mpl_path.LINETO,
        mpl_path.CURVE4,
        mpl_path.CURVE4,
        mpl_path.CURVE4,
        mpl_path.LINETO,
    ]
    
    inner_rad = (start2 - start1) / 2
    outer_rad = (end2 - end1) / 2
    A = (angle2 + angle1) / 2
    L = 4 * np.tan(A/4) / 3
    verts = [
        (x1, y1),
        (x1 + inner_rad * L * np.sin(angle1), y1 - inner_rad * L * np.cos(angle1)),
        (x2 - inner_rad * L * np.sin(angle2), y2 + inner_rad * L * np.cos(angle2)),
        (x2, y2),
        (x3, y3),
        (x3 - outer_rad * L * np.sin(angle2), y3 + outer_rad * L * np.cos(angle2)),
        (x4 + outer_rad * L * np.sin(angle1), y4 - outer_rad * L * np.cos(angle1)),
        (x4, y4),
        (x1, y1)
    ]

    path = mpl_path(verts, codes)
    patch = patches.PathPatch(path, facecolor=facecolor, edgecolor=edgecolor, alpha=0.4)
    ax.add_patch(patch)


def sankey(
    ax, 
    flux_df, 
    angle1=-150,
    angle2=-30,
    row_palette=sns.mpl_palette("tab10", 5),
    col_palette=sns.mpl_palette("tab20", 12),
    sender_length=1,
    receiver_length=1,
    rect_thickness=0.1,
    split_dist=0.,
    pad=0.01,
    add_labels=True,
    add_edges=False,
    norm=True
):
    
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 0.1)
    ax.set_axis_off()
    
    m, n = flux_df.shape
    senders = np.linspace(0, sender_length, m+1)
    receivers = np.linspace(0, receiver_length, n+1)

    for i, j, c, l in zip(senders, senders[1:], row_palette, flux_df.index):
        l = l if add_labels else ""
        rectangular_ideogram(
            ax, i, j, 
            angle=angle1, origin_offset=-split_dist, thickness=-rect_thickness, 
            facecolor=c, label=l, label_ha="right", edgecolor="none" if not add_edges else "0.2"
        )
    for i, j, c, l in zip(receivers, receivers[1:], col_palette, flux_df.columns):
        l = l if add_labels else ""
        rectangular_ideogram(
            ax, i, j, angle=angle2, origin_offset=split_dist, thickness=rect_thickness, 
            facecolor=c, label=l, label_ha="left", edgecolor="none" if not add_edges else "0.2"
        )

    sender_props = receiver_props = flux_df
    if norm:
        sender_props = flux_df / flux_df.sum(axis=1).values[:,None]
        receiver_props = flux_df / flux_df.sum(axis=0)
    
    start_pos = dict(zip(flux_df.index, senders))
    start_pos.update(dict(zip(flux_df.columns, receivers)))
    
    for i, s in enumerate(flux_df.index):
        sender_scale = senders[i+1] - senders[i]
        for j, r in enumerate(flux_df.columns):
            receiver_scale = receivers[j+1] - receivers[j]
            sp = sender_props.loc[s,r]
            rp = receiver_props.loc[s,r]
            if (sp == 0) or (rp == 0):
                continue
            
            chord(
                ax,
                start1=-start_pos[s], 
                end1=-start_pos[s] - sp * sender_scale, 
                start2=start_pos[r], 
                end2=start_pos[r] + rp * receiver_scale,
                angle1=angle1,
                angle2=angle2,
                x_offset=split_dist,
                #y_offset=0,
                facecolor=row_palette[i],
                edgecolor="none" if not add_edges else "0.2"
            )
            start_pos[s] += sp*sender_scale
            start_pos[r] += rp*receiver_scale
        
    return
