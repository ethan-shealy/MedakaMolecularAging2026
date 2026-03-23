#!/bin/bash
#SBATCH --job-name=EMseq-pipeline
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=*****@uga.edu
#SBATCH --out=*****/masterLog.%j



## CONFIG - set parameters for project here -------------------------------------------------------------------

## Set project directory - raw fastq files go in a subdirectory titled raw_reads
PROJECT_DIRECTORY='*****'

## Set sample names here - paired end reads should be in format XXX_1.fq.gz, XXX_2.fq.gz for mates.
SAMPLE_IDS='1,2,3,4......'

## Choose run parameters (TRUE or FALSE) - Deduplication is recommended for whole-genome methylation data, but NOT RRBS

DEDUPLICATION='TRUE'





## PROGRAM (NO EDITS BELOW HERE) ------------------------------------------------------------------------------
echo 'Project directory = ' $PROJECT_DIRECTORY
echo
echo 'Samples: ' $SAMPLE_IDS
echo
echo 'Deduplication? ' $DEDUPLICATION


echo
echo
echo "****************************************************************************"
echo "*                                                                          *"
echo "******************           BEGIN PIPELINE           **********************"
echo "*                                                                          *"
echo "****************************************************************************"
echo
echo


if [ ! -d $PROJECT_DIRECTORY/logs ]
then
    mkdir $PROJECT_DIRECTORY/logs
fi

cd $PROJECT_DIRECTORY/scripts


echo
echo
echo "******************              Genome Index    	      ********************"
echo
echo


if [ ! -d $PROJECT_DIRECTORY/genome/Bisulfite_Genome ]
	then
		sbatch --output=$PROJECT_DIRECTORY/logs/indexLog.%j --chdir=$PROJECT_DIRECTORY --mail-user=$USER@uga.edu genome_index.sh | cut -d ' ' -f4 > indexjob.id
		INDEXID=$(cat indexjob.id)
		echo 'Index job ' $INDEXID ' Started...'
		INDEX='TRUE'
	else	
		echo 'Genome index detected...'
		INDEX='FALSE'
fi
	
	
echo
echo
echo "******************             	Trimming    		  ********************"
echo
echo


if [ ! -d $PROJECT_DIRECTORY/QC ]
then
    mkdir $PROJECT_DIRECTORY/QC
fi

# Create directory for trimmed reads
if [ ! -d $PROJECT_DIRECTORY/trimmed_reads ]
then
    mkdir -p $PROJECT_DIRECTORY/trimmed_reads
fi

if [ $INDEX = 'TRUE' ]
	then
		sbatch --array=$SAMPLE_IDS --output=$PROJECT_DIRECTORY/logs/trimLog.%j --chdir=$PROJECT_DIRECTORY --mail-user=$USER@uga.edu --dependency=afterok:$INDEXID trim_array.sh | cut -d ' ' -f4  > trimjob.id
		TRIMID=$(cat trimjob.id)
		echo 'Trim job ' $TRIMID ' Queued...'
	else
		sbatch --array=$SAMPLE_IDS --output=$PROJECT_DIRECTORY/logs/trimLog.%j --chdir=$PROJECT_DIRECTORY --mail-user=$USER@uga.edu trim_array.sh | cut -d ' ' -f4  > trimjob.id
		TRIMID=$(cat trimjob.id)
		echo 'Trim job ' $TRIMID ' Running...'
fi


echo
echo
echo "******************               Alignment       	    **********************"
echo
echo


#create directory for alignments
if [ ! -d $PROJECT_DIRECTORY/bismark_alignments ]
then
    mkdir -p $PROJECT_DIRECTORY/bismark_alignments
fi


sbatch --array=$SAMPLE_IDS --output=$PROJECT_DIRECTORY/logs/alignLog.%j  --chdir=$PROJECT_DIRECTORY --mail-user=$USER@uga.edu --dependency=afterok:$TRIMID bismark_align_array.sh| cut -d ' ' -f4 > alignjob.id

ALIGNID=$(cat alignjob.id)

echo 'Align ' job $ALIGNID ' Queued...'


echo
echo
echo "******************           Deduplication      	    **********************"
echo
echo


if [ $DEDUPLICATION = 'TRUE' ]
	then
		# Create target directory
		if [ ! -d $PROJECT_DIRECTORY/deduped_alignments ]
			then
				mkdir -p $PROJECT_DIRECTORY/deduped_alignments
		fi
		# Execute dedup script
		sbatch --output $PROJECT_DIRECTORY/logs/dedupLog.%j --chdir=$PROJECT_DIRECTORY --mail-user=$USER@uga.edu --dependency=afterok:$ALIGNID deduplicate.sh | cut -d ' ' -f4 > dedupjob.id
		DEDUPID=$(cat dedupjob.id)
		echo 'Dedup job ' $DEDUPID ' Queued...'
	else
		echo "**** Deduplication skipped (Not recommended for whole-genome data) *****"
fi


echo
echo
echo "******************             Sorting        	    **********************"
echo
echo


#create directory for alignments
if [ ! -d $PROJECT_DIRECTORY/sorted_alignments ]
then
    mkdir -p $PROJECT_DIRECTORY/sorted_alignments
fi

if [ $DEDUPLICATION = 'TRUE' ]
	then
		sbatch --array=$SAMPLE_IDS --output=$PROJECT_DIRECTORY/logs/sortLog.%j  --chdir=$PROJECT_DIRECTORY/deduped_alignments --mail-user=$USER@uga.edu --dependency=afterok:$DEDUPID sort.sh| cut -d ' ' -f4 > sort.id
		SORTID=$(cat sort.id)
		echo 'Sort job ' $SORTID ' Queued...'
	else
		sbatch --array=$SAMPLE_IDS --output=$PROJECT_DIRECTORY/logs/sortLog.%j  --chdir=$PROJECT_DIRECTORY/bismark_alignments --mail-user=$USER@uga.edu --dependency=afterok:$ALIGNID sort.sh| cut -d ' ' -f4 > sort.id
		SORTID=$(cat sort.id)
		echo 'Sort job ' $SORTID ' Queued...'
fi


echo
echo
echo "******************             Methylkit        	    **********************"
echo
echo

if [ ! -d $PROJECT_DIRECTORY/methylkit ]
then
    mkdir -p $PROJECT_DIRECTORY/methylkit
fi


sbatch --array=$SAMPLE_IDS --output=$PROJECT_DIRECTORY/logs/methylkitLog.%j  --chdir=$PROJECT_DIRECTORY --mail-user=$USER@uga.edu --dependency=afterok:$SORTID methylkit.sh


echo
echo
echo "****************************************************************************"
echo "*                                                                          *"
echo "******************           FINISH PIPELINE          **********************"
echo "*                                                                          *"
echo "****************************************************************************"
echo  
echo 