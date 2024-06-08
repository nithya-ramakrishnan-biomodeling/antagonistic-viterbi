#!/bin/bash
#SBATCH --job-name=filter
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=1                                                                              

srun python3 -u filter_binarize.py

