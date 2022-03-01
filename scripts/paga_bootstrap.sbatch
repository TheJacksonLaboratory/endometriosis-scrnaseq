#!/usr/bin/env bash
### SLURM HEADER
#SBATCH --job-name=paga-bootstrap
#SBATCH --output=log/paga-bootstrap-%A-%a.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=bill.flynn@jax.org

#SBATCH --account=robson-lab
#SBATCH --qos=batch
#SBATCH --time=00:25:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB

#SBATCH --export=NONE #ALL
#SBATCH --array=0-99
### SLURM HEADER

group_name=$1
input_data=$2
sample_type=$3

eval "$(conda shell.bash hook)"

localcores=${SLURM_CPUS_ON_NODE}
cd ${SLURM_SUBMIT_DIR}

conda activate endometriosis-2021-12

k=${SLURM_ARRAY_TASK_ID}
outdir=data/bootstrap

python scripts/_paga_bootstrap.py \
    $group_name $input_data $sample_type \
    -k $k \
    -o $outdir \
    -c 50 -n 25 -p 17
