import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjust_text import adjust_text


def volcano_plot(input_list, name, volcano_genes=[], lfc=0.99, pv=1e-2, filename=None):
    dge_list = input_list[~input_list.index.str.startswith("mt-")].copy()
    dge_list["nlogp"] = -np.log10(dge_list.FDR)
    reds = dge_list[(dge_list.FDR < pv) & (dge_list.logFC < -lfc)]
    blues = dge_list[(dge_list.FDR < pv) & (dge_list.logFC > lfc)]
    grays = dge_list[(dge_list.FDR > pv) | (dge_list.logFC.abs() < lfc)]
    volcano_genes = []
    # if len(reds) < 100:
    #    volcano_genes += reds.index.tolist()
    volcano_genes = reds.sort_values("nlogp", ascending=False).index[:30].tolist()
    volcano_genes += reds.sort_values("logFC", ascending=True).index[:4].tolist()

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.grid(False)
    ax.text(0.05, 0.925, "Up in mut", transform=ax.transAxes, fontsize="small")
    ax.text(0.95, 0.925, "Up in wt", ha="right", transform=ax.transAxes, fontsize="small")
    ax.axvline(0, lw=0.5, color="0.4")
    ax.set_title(name)

    scatter_params = dict(edgecolors="none", zorder=10)
    ax.scatter(grays.logFC, grays.nlogp, s=8, c="0.8", alpha=0.5, **scatter_params)
    ax.scatter(blues.logFC, blues.nlogp, s=12, c=sns.xkcd_rgb["medium blue"], alpha=0.9,
               **scatter_params)
    ax.scatter(reds.logFC, reds.nlogp, s=12, c=sns.xkcd_rgb["light red"], alpha=0.9,
               **scatter_params)
    # ax.set_xlim(-5, 5)
    # ax.set_ylim(0, 100)
    ax.set_xlabel("Log Fold Change")
    ax.set_ylabel(r"$-\log_{10}$(q-value)")
    sns.despine(fig, ax)

    labels = []
    for gene in volcano_genes:
        x, y = dge_list.loc[gene, ["logFC", "nlogp"]]
        labels.append(ax.text(x, y, gene, fontsize=5, ha="left", zorder=40))

    if labels:
        adjust_text(
            labels,
            x=reds.logFC.values,
            y=reds.nlogp.values,
            ax=ax, autoalign='xy',
            # force_points=(0.15, 0.15),
            expand_points=(2, 2),
            # expand_text=(1.15, 1.15),
            # force_text=(-0.5, 0.5),
            arrowprops=dict(arrowstyle="-", color='k', lw=0.4)
        )
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=400)
