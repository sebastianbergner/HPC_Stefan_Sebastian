#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name osu_plot
#SBATCH --output=output_plot.log
#SBATCH --exclusive
#SBATCH --time 10

# load modules
module load Anaconda3/2023.10/python-3.11.5-anaconda-2023.10

# plot data
echo "python3 plot.py output_osu.log"
python3 plot.py output_osu.log