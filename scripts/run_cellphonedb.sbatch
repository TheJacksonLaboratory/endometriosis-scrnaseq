#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --job-name=cellphonedb
#SBATCH --output=logs/cellphonedb-custom-%j.log
#SBATCH --mail-type=FAIL

#SBATCH --account=robson-lab
#SBATCH --qos=batch
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=36GB

### SLURM HEADER
usage=$(cat <<EOF
Run with:
  sbatch scripts/run_cellphonedb.sbatch <input_file_prefix>

Expects input data to be in a '\${SLURM_SUBMIT_DIR}/analysis/cellphone' and of the form
  '<prefix>_metadata.txt'
  '<prefix>_counts.txt'
EOF
)

set -euo pipefail

localcores=${SLURM_CPUS_ON_NODE}

cd ${SLURM_SUBMIT_DIR}/analysis/cellphone

file_prefix=${1}
out_dir="${file_prefix}-output"
if [[ -z ${file_prefix} ]]; then
    echo "Need to specify an file_prefix" >&2
    echo $usage >&2
    exit 2
else
    echo "You supplied file prefix: [${file_prefix}]"
    echo "Will store outputs in: [${out_dir}]"
fi

# Load cellphonedb env
module use -a /projects/robson-lab/modules/sumner
module load conda
source activate cpdb

mkdir -p ${out_dir}

cellphonedb method statistical_analysis \
    ${file_prefix}_metadata.txt \
    ${file_prefix}_counts.txt \
    --output-path=${out_dir} \
    --counts-data=gene_name \
    --threads=8 \
    --verbose