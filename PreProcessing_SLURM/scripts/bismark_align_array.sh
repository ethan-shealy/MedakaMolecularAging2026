#!/bin/bash
#SBATCH --job-name=array_align
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --mem=420g
#SBATCH --cpus-per-task=8
#SBATCH --time=4-0
#SBATCH --mail-type=FAIL

## load modules - bismark requires bowtie2
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load Bismark/0.22.3-foss-2019b


echo
echo 'Directory: ' $PWD
echo

## trim files in array style, with values in array index corresponding to sample IDs - e.g. 10_val_1.fastq.gz is sample 10 forward reads, 10_val_2.fastq.gz is reverse, etc.
bismark genome --multicore 8 --maxins 1000 -1 $PWD/trimmed_reads/$SLURM_ARRAY_TASK_ID*val_1.fq.gz -2 $PWD/trimmed_reads/$SLURM_ARRAY_TASK_ID*val_2.fq.gz -o bismark_alignments

## go to directory with newly trimmed/aligned files and run multiqc

ml multiqc/1.11-GCCcore-8.3.0-Python-3.8.2

multiqc -f -n summary_report -o $PWD/QC/ QC bismark_alignments
