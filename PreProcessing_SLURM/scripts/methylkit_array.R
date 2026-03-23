library(methylKit)

# get the input passed from the shell script

array.id <- commandArgs(trailingOnly = TRUE)[1]
print(array.id)

PWD <- commandArgs(trailingOnly = TRUE)[2]
print(PWD)


## Go to directory with sorted and indexed BAM files
setwd(paste0(PWD,"/sorted_alignments"))


# use shell input to paste file name together
file_to_read <- paste0(array.id, "_1_val_1_bismark_bt2_pe.deduplicated_sorted.bam")

print(paste0("Input file:  ", file_to_read))

processBismarkAln(location = file_to_read, sample.id = array.id, assembly = "...", 
                  save.folder = paste0(PWD,"/methylkit"), save.context = c("CpG"), read.context = c("CpG"),
                  mincov = 5, minqual = 20, treatment = NULL)


