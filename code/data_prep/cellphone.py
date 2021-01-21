from pathlib import Path
from typing import Union, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse


def export_for_cellphonedb(
    adata: sc.AnnData,
    outdir: Union[str, Path],
    *,
    partition_key: str = "patient_id",
    partition_prefix: str = "",
    celltype_key: str = "celltype",
    celltype_label: str = "cell_type",
    use_raw: bool = True,
    layer: Optional[str] = None
):
    outdir = Path(outdir)
    if not outdir.exists():
        outdir.mkdir()

    assert (not use_raw) or (layer is None), \
        "Cannot use 'use_raw' and 'layer=' at the same time"

    rownames = adata.var_names
    if use_raw:
        x = adata.raw.X
        rownames = adata.raw.var_names
    elif layer is not None:
        x = adata.layers[layer]
    else:
        x = adata.X

    if issparse(x):
        x = x.todense()

    for partition_id, df in adata.obs.groupby(partition_key):
        metadata = df[celltype_key]
        metadata.index.name = "Cell"
        metadata.columns = [celltype_label]

        inds = adata.obs_names.isin(df.index)
        expr_matrix = pd.DataFrame(
            x[inds, :].T,
            columns=df.index,
            index=rownames
        )

        metadata_path = outdir / f"{partition_prefix}{partition_id}_metadata.txt"
        counts_path = outdir / f"{partition_prefix}{partition_id}_counts.txt"
        output = outdir / f"{partition_prefix}{partition_id}-output"

        metadata.to_csv(metadata_path, sep="\t")
        expr_matrix.to_csv(counts_path, sep="\t")
        output.mkdir(exist_ok=True)


def load_cellphonedb_results(
    adata: sc.AnnData,
    results_dir: Union[str, Path],
    *,
    partition_key: str = "patient_id",
    partition_prefix: str = "",
    celltype_key: str = "celltype",
    celltype_label: str = "cell_type",
    use_raw: bool = True,
    layer: Optional[str] = None
):
    results_dir = Path(results_dir)
    assert results_dir.exists(), f"{results_dir} doesn't exist!"

    partitions = adata.obs[partition_key].unique()

    results = []
    for subdir in results_dir.glob("*-output"):
        if not subdir.is_dir(): continue
        base_name = subdir.stem
        result = curate_ppi_lists(subdir.resolve())
        result["partition_key"] = base_name
        results.append(result)

    return pd.concat(results, axis=0)


def curate_ppi_lists(results_dir):
    results_dir = Path(results_dir)
    sig_means = pd.read_table(results_dir / "significant_means.txt")
    pvals = pd.read_table(results_dir / "pvalues.txt", index_col=1)

    sig_means.drop([
        "id_cp_interaction", "partner_a", "partner_b"
    ], inplace=True, axis=1)

    combos = sig_means.columns[9:]

    rows = []
    for _, srow in sig_means.iterrows():
        meta = srow[:9].values.tolist()
        for combo in combos:

            if not np.isnan(srow[combo]):
                pval = pvals.loc[meta[0], combo]
                row = meta + combo.split("|") + [srow[combo], pval]
                rows.append(row)

    curated = pd.DataFrame(
        rows,
        columns=sig_means.columns[:9].tolist() + ["celltype_a", "celltype_b", "mean",
                                                  "pval"]
    )

    return curated

