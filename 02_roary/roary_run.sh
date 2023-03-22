#!/bin/bash
#
#SBATCH -p lu32
#SBATCH -N 1
#SBATCH --tasks-per-node=8
#SBATCH -t 06:00:00
#SBATCH -J roary_run
#SBATCH -o roary_run_%j.out
#SBATCH -e roary_run_%j.err
#SBATCH --mail-user=yi3446su-s@student.lu.se
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

# write this script to stdout-file useful for scripting errors
cat $0

# load required modules for roary to work
module load GCC/10.3.0 
module load OpenMPI/4.1.1 
module load Roary/3.13.0

# run roary: 8 threads, verbose, use mafft 
roary -p 8 -v -e --mafft *.gff

# clear modules


##################################################################
