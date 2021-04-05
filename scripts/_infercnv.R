#!/usr/bin/env Rscript

library(infercnv)

args = commandArgs(trailingOnly=TRUE)
data_dir <- args[1]
print(args)
print(data_dir)

setwd(data_dir)

counts_file <- paste(data_dir, "counts.txt", sep="/")
gene_file <- paste(data_dir, "gene_pos2.txt", sep="/")
anno_file <- paste(data_dir, "anno.txt", sep="/")
out_dir <- paste(data_dir, "outputs-no-hmm", sep="/")

infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix=counts_file,
    annotations_file=anno_file,
    delim="\t",
    gene_order_file=gene_file,
    ref_group_names=c("normal")
)

infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff=0.1,
    out_dir=out_dir,
    cluster_by_groups=TRUE,
    denoise=TRUE,
    no_plot=FALSE
    HMM=FALSE
)