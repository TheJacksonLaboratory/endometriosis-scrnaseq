#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv


def run_scvelo_paga(_adata, n_pcs=30, n_neighbors=30):
    adata_copy = _adata.copy()
    scv.pp.neighbors(adata_copy, n_pcs=n_pcs, n_neighbors=n_pcs,  metric="correlation")
    scv.pp.moments(adata_copy, n_pcs=n_pcs, n_neighbors=n_neighbors, use_highly_variable=True)
    scv.tl.recover_dynamics(adata_copy)
    top_genes = adata_copy.var['fit_likelihood'].sort_values(ascending=False).index[:100]
    scv.tl.differential_kinetic_test(adata_copy, var_names=top_genes, groupby='subtypes')
    scv.tl.velocity(adata_copy, mode="dynamical", diff_kinetics=True)
    scv.tl.velocity_graph(adata_copy)
    scv.tl.rank_velocity_genes(adata_copy, groupby="subtypes", min_corr=.3)
    scv.tl.paga(adata_copy, groups='subtypes')
    return adata_copy

def subsample_adata(adata, grouping_key="PID", n_cells=50, allow_smallest=True, with_replacement=True):
    min_group_size = adata.obs[grouping_key].value_counts().min()
    
    if allow_smallest:
        if min_group_size < n_cells:
            print(f"Warn: smallest group has {min_group_size} which is < {n_cells}.  Using that number instead.")
        n_cells = min(n_cells, min_group_size)

    idx = []
    for group in adata.obs.groupby(grouping_key).groups.values():
        if len(group) < (n_cells / 4): continue
        idx += np.random.choice(group, size=n_cells, replace=with_replacement).tolist()
    return adata[idx, :].copy()

def bootstrap(_adata, sample_type, n_cells, n_neighbors, n_pcs):
    _boot = subsample_adata(
        _adata[_adata.obs.sample_type_rename.isin(sample_type)], 
        grouping_key="subtypes", n_cells=n_cells, allow_smallest=False 
    )
    del _boot.uns["neighbors"]
    del _boot.obsp
    boot = run_scvelo_paga(_boot, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return boot.uns["paga"]


def main(args):
    adata = sc.read(args.h5ad_path)

    if not args.outdir.exists():
        args.outdir.mkdir(parents=True)

    paga_dict = bootstrap(
        adata, args.sample_type, args.n_cells, args.n_neighbors, args.n_pcs
    )

    prefix = f"{args.name}_{'-'.join(args.sample_type)}"
    suffix = f"c{args.n_cells}_n{args.n_neighbors}_p{args.n_pcs}_k{args.k}.txt"
    fname = "_".join((prefix, "{}", suffix))

    keys = ["connectivities", "connectivities_tree", "transitions_confidence"]
    for key in keys:
        np.savetxt(args.outdir / fname.format(key), paga_dict[key].todense())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("name", type=str)
    parser.add_argument("h5ad_path", type=Path)
    parser.add_argument("sample_type", type=str, nargs="+")
    parser.add_argument("-k", "--k", type=int)
    parser.add_argument("-c", "--n-cells", type=int, default=50)
    parser.add_argument("-n", "--n-neighbors", type=int, default=25)
    parser.add_argument("-p", "--n-pcs", type=int, default=17)
    parser.add_argument("-o", "--outdir", type=Path)
    args = parser.parse_args()
    main(args)
