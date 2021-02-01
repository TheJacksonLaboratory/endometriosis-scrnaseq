import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjust_text import adjust_text


def volcano_plot(
    input_list, title, lfc=0.99, pv=1e-2, neg_label="", pos_label="",
    label_genes=False, genes_to_label=[], top_neg_genes=10, top_pos_genes=10,
    filename=None, adjust_kws={}
):
    dge_list = input_list[~input_list.index.str.startswith("mt-")].copy()
    dge_list["nlogp"] = -np.log10(dge_list.FDR)
    reds = dge_list[(dge_list.FDR < pv) & (dge_list.logFC < -lfc)]
    blues = dge_list[(dge_list.FDR < pv) & (dge_list.logFC > lfc)]
    grays = dge_list[(dge_list.FDR > pv) | (dge_list.logFC.abs() < lfc)]

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.grid(False)
    ax.text(0.05, 0.925, neg_label, transform=ax.transAxes, fontsize="small")
    ax.text(0.95, 0.925, pos_label, ha="right", transform=ax.transAxes, fontsize="small")
    ax.axvline(0, lw=0.5, color="0.4")
    ax.set_title(title)

    scatter_params = dict(edgecolors="none", zorder=10)
    ax.scatter(grays.logFC, grays.nlogp, s=8, c="0.8", alpha=0.5, **scatter_params)
    ax.scatter(blues.logFC, blues.nlogp, s=12, c=sns.xkcd_rgb["medium blue"], alpha=0.9,
               **scatter_params)
    ax.scatter(reds.logFC, reds.nlogp, s=12, c=sns.xkcd_rgb["light red"], alpha=0.9,
               **scatter_params)
    ax.set_xlabel("Log Fold Change")
    ax.set_ylabel(r"$-\log_{10}$(q-value)")
    sns.despine(fig, ax)

    if label_genes:
        for (df, n_top, asc) in zip(
            (reds, blues),
            (top_neg_genes, top_pos_genes),
            (True, False)
        ):
            gene_labels = []
            if len(genes_to_label) == 0:
                genes = df.sort_values(by="logFC", ascending=asc).index[:n_top].tolist()
                for gene in genes:
                    x, y = dge_list.loc[gene, ["logFC", "nlogp"]]
                    gene_labels.append(ax.text(x, y, gene, fontsize=5, ha="left", zorder=40))
            else:
                for gene in genes_to_label:
                    if gene not in df.index: continue
                    x, y = df.loc[gene, ["logFC", "nlogp"]]
                    gene_labels.append(ax.text(x, y, gene, fontsize=5, ha="left", zorder=40))

            default_kws = dict(
                autoalign='xy',
                expand_points=(2, 2),
                expand_text=(1.15, 1.15),
                arrowprops=dict(arrowstyle="-", color='k', lw=0.4)
            )
            default_kws.update(adjust_kws)
            adjust_text(
                gene_labels,
                x=df.logFC.values,
                y=df.nlogp.values,
                ax=ax, **default_kws
            )

    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=400)


if __name__ == "__main__":
    print('hi')