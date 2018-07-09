#!/bin/bash

#SBATCH --job-name=execute_RNAseq_pipeline 
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=0:29:00   # modify this to reflect how long to let the job go. This indicates 4 hours.
#SBATCH --output=log_RNAseq_pipe_%j.txt

##
source /projects/dcking@colostate.edu/paths.bashrc

## execute the RNA-seq_pipeline
#bash RNAseq_analyzer_mouse_180706.sh ../04_testing/metadata_mouse.txt 12
     # modify the SECOND argument to point to YOUR metadata.file
     # modify the THIRD argument to indicate the number of THREADS you 
     # want to use. This number must match the number in #SBATCH --ntasks=#

## clean up by zipping .fastq files and deleting extra files
bash RNAseq_cleanup_mouse_180706.sh ../04_testing/metadata_mouse.txt
     # modify the SECOND argument to point to YOUR metadata.file
