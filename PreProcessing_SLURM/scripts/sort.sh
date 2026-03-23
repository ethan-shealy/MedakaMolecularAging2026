#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=3-0
#SBATCH --mail-type=FAIL


module load SAMtools/1.16.1-GCC-11.3.0


echo
echo 'Directory: ' $PWD
echo


for i in ./*.bam
	do  samtools sort -O BAM -@ 4 $i > ../sorted_alignments/${i/.bam/_sorted.bam} 
done

