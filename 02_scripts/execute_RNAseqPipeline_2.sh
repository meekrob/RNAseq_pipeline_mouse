#!/bin/bash

#SBATCH --job-name=execute_RNAseq_pipeline2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=0:29:00
#SBATCH --output=log_RNAseq_pipe2_%j.txt


##
source /projects/dcking@colostate.edu/paths.bashrc

## execute the RNA-seq_pipeline
#bash RNAseq_analyzer_mouse_180603.sh ../01_input/metadata_mouse.txt
bash RNAseq_cleanup_mouse_180706.sh ../01_input/metadata_mouse.txt
