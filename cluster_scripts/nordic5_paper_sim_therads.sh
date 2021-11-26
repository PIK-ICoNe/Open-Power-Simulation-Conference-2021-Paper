#!/bin/bash

#SBATCH --qos=medium
#SBATCH --partition=standard
#SBATCH --account=icone
#SBATCH --error=%x-%j-%N.err
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --tasks-per-node=20
#SBATCH --job-name=nordic5_multi_100_paper_sim
#SBATCH --mem=64GB

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.6.1
julia --threads=$SLURM_NTASKS nordic5_paper_sim_threads.jl