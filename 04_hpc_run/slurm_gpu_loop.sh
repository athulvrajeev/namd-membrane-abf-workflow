#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=40
#SBATCH --time=96:00:00
#SBATCH --job-name=ABF_21-40
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --partition=gpu
#SBATCH --gres=gpu:2

cd $SLURM_SUBMIT_DIR

module purge
module load apps/namd/2.14/gpu

# 2 MPI ranks per node
NAMDMPI=$(expr $SLURM_NTASKS / 20)

for ((i=20; i<=39; i++)); do
    folder="win$i"
    conf="${folder}/${folder}.conf"  # Configuration file inside the folder
    output="${folder}/${folder}.log" # Log file inside the folder
    
    echo "Running simulation for $folder"
    (time mpiexec.hydra -n $NAMDMPI namd2 +ppn 19 +pemap 1-19,21-39 +commap 0,20 +setcpuaffinity +isomalloc_sync +idlepoll $conf +devices 0,1) >& $output
    wait
    echo "Simulation for $folder finished."
done

