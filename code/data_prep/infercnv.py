from typing import Union
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse


def save_for_infercnv(
    adata: sc.AnnData,
    normal_vs_other_key: str = "sample_status",
    normal_value: str = "normal",
    grouping_key: str = "patient_id",
    layer: Union[str, None] = None,
    output_dir: Union[Path, str, None] = None
):
    """
    This function takes an adata object and writes the counts and annotation file
    needed by inferCNV.
    :param adata: AnnData object
    :param normal_vs_other_key: .obs key to distinguish normal cells from non-normal cells
    :param normal_value: value in `normal_vs_other_key` that denotes normal cells
    :param grouping_key: .obs key that is to be used to add additional metadata to non-normal cells
    :param layer: the layer to use for expression. if None, will try to use adata.raw.X if it exists otherwise adata.X
    :param output_dir: the directory under which data will be saved. Must exist prior to execution
    """

    annos = pd.DataFrame(index=adata.obs_names)
    ref_values = adata.obs[normal_vs_other_key]
    annos["label"] = ref_values.map({normal_value: "normal"}).fillna(ref_values)
    annos["label"] += "_" + adata.obs[grouping_key]

    if output_dir is None:
        return annos
    else:
        output_dir = Path(output_dir)
        assert output_dir.exists(), f"Dir {output_dir} doesn't exist!"

    X = adata.X
    if layer is not None:
        X = adata.layers[layer]
    elif adata.raw is not None:
        X = adata.raw.X

    X = X.T.todense() if issparse(X) else X.T

    counts_matrix = pd.DataFrame(
        X,
        index=adata.var["gene_ids"],
        columns=adata.obs_names
    )
    counts_matrix.to_csv(output_dir / "counts.txt", sep="\t")

    annos.to_csv(
        output_dir / "anno.txt", sep="\t", header=False
    )


def test_data():
    adata = sc.AnnData(
        np.arange(60).reshape(10, 6),
        var=pd.DataFrame(index=[f"gene_{k}" for k in "abcdef"]),
        obs=pd.DataFrame(index=[f"cell_{k}" for k in range(1, 11)]),
    )
    adata.obs["norm_v_tumor"] = ["norm"]*3 + ["tumor-A"]*2 + ["tumor-B"]*5
    adata.obs["patient_id"] = [f"P_{k}" for k in np.random.randint(1, 3, size=10)]
    return adata


if __name__ == "__main__":
    save_for_infercnv(
        test_data(),
        "norm_v_tumor",
        "norm",
        "patient_id"
    )
