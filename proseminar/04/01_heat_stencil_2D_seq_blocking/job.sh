#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps4.optim
#SBATCH --output=output/output.log
#SBATCH --open-mode=append
#SBATCH --time=00:10:00
#SBATCH --mail-user=sebastian.bergner@student.uibk.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
# stefan.r.wagner@student.uibk.ac.at,
#SBATCH --exclusive

echo "$@"
eval "$@"
echo
