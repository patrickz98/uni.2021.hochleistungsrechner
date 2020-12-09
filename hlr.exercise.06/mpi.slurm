#!/bin/bash

# Run a total of three tasks on two nodes (that is, five tasks per node).
#SBATCH --nodes=2
#SBATCH --ntasks=3
#SBATCH --partition=west
#SBATCH --output=job.out
#SBATCH --error=job.err

mpiexec ./timempi
