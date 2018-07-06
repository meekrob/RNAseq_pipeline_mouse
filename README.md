# RNAseq_pipeline_mouse
This is a very basic RNA-seq pipeline that I use for analyzing mouse RNA-seq paired-end sequencing. Step1 is a simple wrapper that performs quality control, genome alignment, basic format conversions, and htseq-count tabulation for paired-end RNA-seq samples using the mouse genome. Step2 is a clean up program that removes unnecessary files and compress files to save space.

### Programs:
execute_RNAseq_pipeline.sh
RNAseq_analyzer_mouse_180706.sh
RNAseq_cleanup_mouse_180706.sh

### Author:
Erin Osborne Nishimura

### Date initiated:
July 3, 2018

### Dependencies: 
  Requires fastqc, hisat2, htseq, samtools, deep=tools

### Requires: 
1) INPUT: .fastq.gz files. For each sample, paired forward and reverse sequencing files are required. These should be placed in an input directory.
2) INPUT: \_metadata.txt file: A metadata file with two columns. The first two columns are fastq.gz file names. The third column is a "nickname" of each sample. Later columns can be included with other metadata information. Metadata file should be placed within the inputdir directory. Example of a metadata file:
3) BUILD: 
#
# Requires: bt2 files for the mouse genome.
# 
# Requires a gtf file for the mouse genome.
#
# Requires a .fa genome sequence for the mouse genome.
# 
# Executed with:
# bash RNAseq_analyzer_mouse.sh metadata.txt 2>&1 | tee 180423_output.txt
################################################
