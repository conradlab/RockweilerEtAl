#!/bin/bash

#SBATCH -o realigned_bam2filtered_readcount.o%J
#SBATCH -e realigned_bam2filtered_readcount.e%J
#SBATCH --mem=8G
# DESCRIPTION: slurm wrapper for realigned_bam2filtered_readcount.bash 

bam=$1
output_dir=$2

bash /home/nrockweiler/twins_uk/src/alignment/realigned_bam2filtered_readcount.COPY.COPY.bash $bam $output_dir
