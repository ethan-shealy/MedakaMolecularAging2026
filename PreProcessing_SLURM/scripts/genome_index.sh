#!/bin/bash
#SBATCH --job-name=index
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL


module load Bowtie2/2.3.5.1-GCC-8.3.0
module load Bismark/0.22.3-foss-2019b

bismark_genome_preparation genome
