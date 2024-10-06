#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name osu_plot
#SBATCH --output=logs/plot.log
#SBATCH --exclusive
#SBATCH --time 10

# load modules
module load Anaconda3/2023.10/python-3.11.5-anaconda-2023.10

# plot data
echo "jobs/plot.py logs/01_osu_different_cores_of_same_socket.log logs/02_osu_different_sockets_of_same_node.log logs/03_osu_different_nodes.log"
python3 jobs/plot.py logs/01_osu_different_cores_of_same_socket.log logs/02_osu_different_sockets_of_same_node.log logs/03_osu_different_nodes.log