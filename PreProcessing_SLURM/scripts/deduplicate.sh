#!/bin/bash
#SBATCH --job-name=dedup
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=3-0
#SBATCH --mail-type=FAIL

module load Bismark/0.22.3-foss-2019b
module load SAMtools/1.16.1-GCC-11.3.0


echo
echo 'Directory: ' $PWD
echo

deduplicate_bismark -p --bam --output_dir $PWD/deduped_alignments $PWD/bismark_alignments/*.bam