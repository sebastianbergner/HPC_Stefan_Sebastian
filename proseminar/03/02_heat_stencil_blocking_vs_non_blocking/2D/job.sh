#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps3.2_2D
#SBATCH --output=output/output.log
#SBATCH --open-mode=append
#SBATCH --time=00:10:00
#SBATCH --exclusive

echo "$@"
eval "$@"
echo
