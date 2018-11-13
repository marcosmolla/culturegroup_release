#!/bin/bash
#SBATCH --job-name=learn
#SBATCH --array=1-30%5

srun julia run.jl
