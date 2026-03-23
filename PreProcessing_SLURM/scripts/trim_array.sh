#!/bin/bash
#SBATCH --job-name=array_trim
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G
#SBATCH --time=3-0
#SBATCH --mail-type=FAIL


##Load modules (pigz and Python 3 required for multicore TrimGalore)
ml Python/3.9.5-GCCcore-10.3.0
ml pigz/2.6-GCCcore-10.2.0
ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

echo
echo 'Directory: ' $PWD
echo

## trim files in array style, with values in array index corresponding to sample IDs - e.g. input10_1.fastq.gz is sample 10 forward reads, input10_2.fastq.gz is reverse, etc.
trim_galore $PWD/raw_reads/${SLURM_ARRAY_TASK_ID}_{1,2}.fq.gz --cores 2 --quality 20 --paired --illumina --clip_r1 8 --clip_r2 8 --three_prime_clip_r1 8 --three_prime_clip_r2 8 --fastqc --fastqc_args "--outdir QC" -o trimmed_reads
