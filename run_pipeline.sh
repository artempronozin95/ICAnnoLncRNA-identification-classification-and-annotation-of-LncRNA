#!/bin/bash

#SBATCH -A pronozinau
#SBATCH --mem 100GB
#SBATCH -p common
#SBATCH -N 1
#SBATCH -t 144:00:00

snakemake -j 1 --use-conda