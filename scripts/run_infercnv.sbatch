#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --job-name=infercnv
#SBATCH --output=/projects/robson-lab/research/endometriosis/logs/infercnv-%j.log
#SBATCH --mail-type=FAIL

#SBATCH --qos=batch
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=96000

#SBATCH --export=NONE #ALL
### SLURM HEADER

localcores=${SLURM_TASKS_PER_NODE}

cd ${SLURM_SUBMIT_DIR}

module load singularity
img="/projects/robson-lab/research/endometriosis/tools/jupyterlab-with-r.sif"

data_dir=$1
if [ -z ${data_dir} ]; then
    echo "No directory with infercnv data specified" >&2
    exit 1
elif [ ! -d ${data_dir} ]; then
    echo "Infercnv data dir doesn't exist: ${data_dir}" >&2
    exit 2
fi

#script_dir=$(realpath $(dirname $0))
script_dir="/projects/robson-lab/research/endometriosis/manuscript/scripts"
subscript="${script_dir}/_infercnv.R"

# HACKKKK
cp /projects/robson-lab/research/ibc/analysis/infercnv/IBC-epithelial/inputs/gene_pos2.txt ${data_dir}

singularity exec ${img} Rscript ${subscript} ${datadir}