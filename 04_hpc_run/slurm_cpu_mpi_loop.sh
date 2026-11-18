#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks-per-node=48
#SBATCH --time=96:00:00
#SBATCH --job-name=athul_md
#SBATCH --error=job.%J.err_2_node_48
#SBATCH --output=job.%J.out_2_node_48
#SBATCH --partition=cpu

# Load the Open MPI module
module load openmpi

# Change directory to your working directory
cd $SLURM_SUBMIT_DIR
source ~/.bashrc  # <-- This ensures namd2 is in your PATH

# Run simulations for win0 to win20 sequentially
for ((i=0; i<=19; i++)); do
    folder="win$i"
    conf="${folder}/${folder}.conf"  # Configuration file inside the folder
    output="${folder}/${folder}.log" # Log file inside the folder
    
    echo "Running simulation for $folder"
    time mpirun -n $SLURM_NTASKS  namd2 $conf > $output
    wait
    echo "Simulation for $folder finished."
done

