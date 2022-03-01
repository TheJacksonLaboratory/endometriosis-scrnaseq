#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from pathlib import Path
from datetime import datetime


def save_adata(adata, suffix="", subdir=""):
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d")

    filename = f"{adata.uns['sampleid']}{'-' + suffix if suffix else ''}-{timestamp}.h5ad"
    sc.write(Path(adata.uns["output_dir"]) / subdir / filename, adata)


def save_for_loupe(adata, outdir, cluster_key="cluster_revised"):
    final_barcodes = adata.obs_names
    loupe_barcodes = final_barcodes.str.extract("([ACGT]+)", expand=True)
    loupe_barcodes["num"] = final_barcodes.str.extract("[ACGT]+-\d+-(\d+)", expand=False)
    loupe_barcodes["num"] = (loupe_barcodes["num"].astype(int) + 1).astype(str)
    loupe_barcodes = pd.DataFrame(loupe_barcodes.iloc[:,0] + "-" + loupe_barcodes.iloc[:,1], columns=["Barcode"])
    cloupe_umap = pd.DataFrame(
        adata.obsm["X_umap"],
        index=loupe_barcodes.Barcode.values,
        columns=["x coordinate", "y coordinate"]
    )
    
    cloupe_umap.index.name = "barcode"
    print(cloupe_umap.head())

    cloupe_metadata = pd.DataFrame({
        "Custom Cluster": "Cluster " + adata.obs[cluster_key].str.zfill(2).values,
        "Celltype": adata.obs["celltype"].values,
        "Patient ID": adata.obs["Patient_id"].values,
        "Fragment ID": adata.obs["sampleid"].values,
        "Sample Type":adata.obs["sample_type"].values,
        "Subtypes":adata.obs["subtypes"].values,
        "Old subtypes":adata.obs["old_subtypes"].values
    },
        index=loupe_barcodes.Barcode.values
    )

    cloupe_metadata.index.name = "barcode"
    print(cloupe_metadata.head())
    
    outpath = Path(f'{outdir}/cloupe/{sub_dir}/')
    loupe_barcodes.to_csv(outpath / "lymphocytes-loupe_barcodes.csv", index=False)
    cloupe_metadata.to_csv(outpath / "lymphocytes-loupe_metadata.csv")
    cloupe_umap.to_csv(outpath / "lymphocytes-loupe_umap.csv")
