# RNAseq_pipeline_mouse
This is a very basic RNA-seq pipeline that I use for analyzing mouse RNA-seq paired-end sequencing. Step1 is a simple wrapper that performs quality control, genome alignment, basic format conversions, and htseq-count tabulation for paired-end RNA-seq samples using the mouse genome. Step2 is a clean up program that removes unnecessary files and compress files to save space.

### Programs:
execute_RNAseq_pipeline.sh
RNAseq_analyzer_mouse_180706.sh
RNAseq_cleanup_mouse_180706.sh

### Author:
Erin Osborne Nishimura

### Date initiated:
July 6, 2018

### Dependencies: 
  Requires fastqc, hisat2, htseq, samtools, deep=tools

### Requires: 
1. INPUT: .fastq.gz files. For each sample, paired forward and reverse sequencing files are required. These should be placed in an input directory.
2. INPUT: \_metadata.txt file: A metadata file with two columns. The first two columns are fastq.gz file names. The third column is a "nickname" of each sample. Later columns can be included with other metadata information. Metadata file should be placed within the inputdir directory. Example of a metadata file:
3. BUILD: .bt2 files for the mouse genome. These are produced using hisat2-build. For instructions see https://ccb.jhu.edu/software/hisat2/manual.shtml#the-hisat2-build-indexer
4. GENOME: .fa file for the mouse genome. This is the sequence of the mouse genome.
5. GENOME: .gtf file for the mouse genome. This is a genome annotation file of gene features. Version and coordinates must match the genome sequence (.fa above)

### Executed on SUMMIT with:
```
$ sbatch execute_RNAseq_pipeline.sh
```

##################

# INSTRUCTIONS

