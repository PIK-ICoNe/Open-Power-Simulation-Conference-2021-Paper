#!/bin/bash

#SBATCH --qos=medium
#SBATCH --partition=standard
#SBATCH --account=Icone
#SBATCH --error=%x-%j-%N.err
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --tasks-per-node=1
#SBATCH --job-name=nordic5_paper_sim
#SBATCH --mem=15GB

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.6.1
julia nordic5_paper_sim.jl $SLURM_NTASKS
