#!/bin/bash

#SBATCH -A pronozinau
#SBATCH --mem-per-cpu=10GB
#SBATCH -p common
#SBATCH -N 1
#SBATCH -t 03:00:00

#rm -R TMHMM_3164369/

python scripts/tissue.py