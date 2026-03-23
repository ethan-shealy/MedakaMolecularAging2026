#!/bin/bash     
#SBATCH --job-name=methylkit_array 
#SBATCH --partition=batch           		                    
#SBATCH --ntasks=1
#SBATCH --mem=100G               
#SBATCH --cpus-per-task=2
#SBATCH --time=1-0
#SBATCH --mail-type=END,FAIL

echo
echo 'Directory: ' $PWD
echo

# Load R
module load R/4.2.1-foss-2020b

### Initialization
# Get Array ID
i=${SLURM_ARRAY_TASK_ID}


# Pass line #i to a R script 
Rscript --vanilla $PWD/scripts/methylkit_array.R ${i} ${PWD}


echo
echo '******************** FINISHED ***********************'
echo
