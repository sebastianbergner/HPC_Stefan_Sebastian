#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps3.1
#SBATCH --output=output/output.log
#SBATCH --open-mode=append
#SBATCH --time=00:10:00
# #SBATCH --mail-user=stefan.r.wagner@student.uibk.ac.at,sebastian.bergner@student.uibk.ac.at
# #SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --exclusive

echo "$@"
eval "$@"
echo
