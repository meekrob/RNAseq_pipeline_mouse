# RNAseq_pipeline_mouse
This is a very basic RNA-seq pipeline that I use for analyzing mouse, paired-end RNA-seq. Step1 is a simple wrapper that performs quality control, genome alignment, basic format conversions, and htseq-count tabulation for paired-end RNA-seq samples using the mouse genome. Step2 is a clean up program that removes unnecessary files and compress files to save space.

### Programs:
Wrapper: execute_RNAseq_pipeline.sh  
Step1: RNAseq_analyzer_mouse_180706.sh  
Step2: RNAseq_cleanup_mouse_180706.sh  

### Author:
Erin Osborne Nishimura

### Date initiated:
July 6, 2018

### Dependencies: 
  Requires fastqc, hisat2, htseq, samtools, deep-tools

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

-----

# INSTRUCTIONS

### Step1: Modify RNAseq_analyzer_mouse_180706.sh
*  Open RNAseq_analyzer_180706.sh and modify the section between ###### MODIFY THIS ####### and ####### DONE MODIFYING ######
* Update __input directory__
  * Modify the this section to point to your input directory where the ```_metadata.txt``` file and ```.fastq.gz``` files can be located.
  ```
  #The input samples (metadata file and _fastq.gz files) live in directory:
  inputdir="../01_input/"
  ```
* Update __.bt2 directory__
  * Modify this to point to the path where the hisat2 generated ```.bt2``` files can be located. The default points to Erin's mm10 .bt2 files.
  ```
  #This is where the bt2 files live:
  hisat2path="/projects/erinnish@colostate.edu/genomes/mm10/from_ucsc/mm10"
  ```
* Update the path to the __genome sequence__
  * Modify this to point to the file of the full mouse genome. The default points to Erin's mm10 version GRCm38/mm10 downloaded from UCSC http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz on 180424.
chromFa.tar.gz          09-Feb-2012 13:54  830M 
  ```
  #This is where the genome sequence lives:
  genomefa="/projects/erinnish\@colostate.edu/genomes/mm10/from_ucsc/chromFa.tar.gz"
  ```
* Update the path to the __genome feature file__
  * Modify this to point to the file listing all the mouse genome features. The default points to Erin's mm10 .gtf file. This was downloaded from... ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/. The original file was Mus_musculus.GRCm38.92.chr.gtf.gz. The file was converted to match the UCSC mm10 genome using a home-made script called ```180420_convert_chromosome_names.awk```. This converted the names of each chromosome as outlined in the Ensembl database to those in the UCSC database. The resulting file is called ```Mus_musculus_GRCm38_2UCSC.gtf```.
  ```
  #This is where the gtf file lives:
  gtffile="/projects/erinnish@colostate.edu/genomes/mm10/from_ensembl/gtf/Mus_musculus_GRCm38_2UCSC.gtf"
  ```
* Specify whether __ERCC Spike-In Controls__ were used 
  * Update the toggle to determine whether ERCC spike-in controls were used to __TRUE__ or __FALSE__
  ```
  ercc="FALSE"   # Change to TRUE if ERCC spike-ins were used in the experiment
  ```
* IF ERCC Spike-ins were used, update the path to the directory containing __.bt2 files for ERCC Spike-in controls__
  ```
  #This is where the ercc bt2 files lives:
  erccpath="/projects/erinnish@colostate.edu/genomes/ercc/ercc92"
  ```
* IF ERCC Spike-ins were used, update the path the __.gtf files for ERCC Spike-in controls__
  ```
  #This is where the ercc .gtf file lives:
  erccgtf="/projects/erinnish@colostate.edu/genomes/ercc/ERCC92.gtf"
  ```



####### MODIFY THIS SECTION #############

#The input samples (metadata file and _fastq.gz files) live in directory:
inputdir="../01_input/"

#This is where the bt2 files live:
hisat2path="/projects/erinnish@colostate.edu/genomes/mm10/from_ucsc/mm10"

#This is where the gtf file lives:
gtffile="/projects/erinnish@colostate.edu/genomes/mm10/from_ensembl/gtf/Mus_musculus_GRCm38_2UCSC.gtf"

#This is where the genome sequence lives:
genomefa="/projects/erinnish\@colostate.edu/genomes/mm10/from_ucsc/chromFa.tar.gz"
         
#Number of threads to use:
pthread=11



# WERE ERCC SPIKE INS USED?
ercc="FALSE"   # Change to TRUE if ERCC spike-ins were used in the experiment

#This is where the ercc bt2 files lives:
erccpath="/projects/erinnish@colostate.edu/genomes/ercc/ercc92"

#This is where the ercc .gtf file lives:
erccgtf="/projects/erinnish@colostate.edu/genomes/ercc/ERCC92.gtf"


########## DONE MODIFYING ###############
