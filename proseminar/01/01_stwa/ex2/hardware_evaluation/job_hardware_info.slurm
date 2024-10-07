#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name lstopo
#SBATCH --output=output_hardware_info.log
#SBATCH --exclusive
#SBATCH --time 10

# load modules
module load hwloc

# run programs
echo "lstopo --of txt"
lstopo --of txt