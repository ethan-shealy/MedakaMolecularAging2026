# EMseq-pipeline
WIP pipeline for the analysis of whole-genome enzymatic methylation sequencing.  


To begin:


1) Place your reference genome in FASTA format into the folder named "Genome". Place gzipped fastq files containing your raw read data into the folder "raw_reads". 
Filenames should be formatted as follows for paired-end data, where XXX is replaced by the numeric identifier for each sample: 
XXX_1.fq.gz for forward reads, XXX_2.fq.gz for reverse reads.


2) Open master.sh - edit the header to send it to your own email, and replace the --output destination with your own directory. 
Under the --array option, list your sample IDs seperated by commas with no spaces (e.g. --array=01,02,03,04,05). 
These numbers should correspond to the filenames you placed in "raw_reads".


3) Move down to the 'Config' section of master.sh and specify where your directory is located using the PROJECT_DIRECTORY='/your/file/location/EMseq-Pipeline' 
(make sure the single quotes are present). Also make sure to fix the SLURM header and specify your email address and directory. Change DEDUPLICATE to 'FALSE' if you
are working with RRBS data. 


4) That's about it - now run master.sh using sbatch. QC files should end up in a QC folder, and log files for each sample will end up in 'logs'. 
The resulting output files will end up in correspondingly named folders as well. 


