#!/usr/bin/env bash

#SBATCH --job-name=execute_RNAseq_pipeline 
#SBATCH --nodes=1
#SBATCH --ntasks=12      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=3:00:00   # modify this to reflect how long to let the job go. This indicates 4 hours.
#SBATCH --output=log_RNAseq_pipe_%j.txt

## Container setup
module load singularity 
CONTAINER_IMG=/projects/dcking@colostate.edu/shub/monaghaa-dcking-master-latest.simg

# Your metadata.file
METADATA=/projects/jesshill@colostate.edu/erins/RNAseq_pipeline_mouse/01_input/metadata_mouse.txt

if true
then
    ## execute the RNA-seq_pipeline
    singularity exec $CONTAINER_IMG \
    bash RNAseq_analyzer_mouse_180706.sh $METADATA $SLURM_NTASKS
else
    ## clean up by zipping .fastq files and deleting extra files
    # does samming and bamming
    singularity exec $CONTAINER_IMG \
    bash RNAseq_cleanup_mouse_180706.sh $METADATA
fi
