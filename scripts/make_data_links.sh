#!/usr/bin/env bash 
set -euo pipefail

cur_file=$(readlink -f $0)
cur_dir="$(dirname "$cur_file")"
par_dir="$(dirname "$cur_dir")"
data_dir="${par_dir}/data"

src_dir="/projects/robson-lab/research/endometriosis/data/h5ad/"
batch_name="Endo-Tissue-EC19001-SC2100428-3"

objs=(
    "${batch_name}-main-final-cc-20220105.h5ad"
    "Myeloid/${batch_name}-mye-final-cc-20220105.h5ad"
    "TNK/${batch_name}-lym-final-cc-20220105.h5ad"
    "Epithelial/${batch_name}-epi-final-cc-20220105.h5ad"
    "Stromal/${batch_name}-stro-final-cc-20220105.h5ad"
    "Endothelial/${batch_name}-endo-final-cc-20220105.h5ad"
    "non_immune/${batch_name}-org-final-cc-20220105.h5ad"
)
dests=(main mye lym epi stro endo org)
for k in $(seq 0 $((${#objs[@]}-1)))
do
    src="${src_dir}/${batch_name}/${objs[$k]}"
    dst="${data_dir}/${dests[$k]}.h5ad"
    ln -s "$src" "$dst"
done
