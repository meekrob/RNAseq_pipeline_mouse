# RNAseq_pipeline_mouse
This is a very basic RNA-seq pipeline that I use for analyzing mouse, paired-end RNA-seq. **Step1** is a simple wrapper that performs quality control, genome alignment, basic format conversions, and htseq-count tabulation for paired-end RNA-seq samples using the mouse genome. **Step2** is a clean up program that removes unnecessary files and compress files to save space.

### Programs:
**Wrapper:** execute_RNAseq_pipeline.sh  
**Step1:** RNAseq_analyzer_mouse_180706.sh  
**Step2:** RNAseq_cleanup_mouse_180706.sh  

### Author:
Erin Osborne Nishimura

### Date initiated:
July 6, 2018

### Dependencies: 
* Requires the installation of the follwing software: 
  * fastqc
  * hisat2
  * htseq
  * samtools
  * deep-tools
* Requires group access to the Nishimura lab on SUMMIT

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

## 1: Modify Step1... RNAseq_analyzer_mouse_180706.sh
*  Navigate into ```02_scripts```
*  Open ```RNAseq_analyzer_180706.sh``` and modify the section between ###### MODIFY THIS ####### **and** ####### DONE MODIFYING ######
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
## 2: Modify Step2... RNAseq_cleanup_mouse_180706.sh
*  Navigate into ```02_scripts```
* Open ```RNAseq_cleanup_180706.sh``` and modify the section between ###### MODIFY THIS ####### **and** ####### DONE MODIFYING ######
* Update __input directory__
  * Modify the this section to point to your input directory where the ```_metadata.txt``` file and ```.fastq.gz``` files can be located.
  ```
  #The input samples (metadata file and _fastq.gz files) live should in directory. Optionally, if you plan on running a testing run first, set this option to "../04_testing/"
  inputdir="../01_input/"
  ```
* Update the **output directory date**
  * If the original analysis was done on a different date, update that date. 
  ```
  #This is the output_directory:
  #DATE=`date +%Y-%m-%d`
  #OR
  DATE=2018-07-06
  outputdir="../03_output/"$DATE"_output/
  ```
  If it was done today, you can opt use system date to auto-input the date, like so...
  ```
  #This is the output_directory:
  DATE=`date +%Y-%m-%d`
  #OR
  #DATE=2018-07-06
  outputdir="../03_output/"$DATE"_output/
  ```
 


 ## 3: Optional testing... execute_RNAseq_pipeline.sh to run on test files
 * Open **execute_RNAseq_pipeline.sh** and modify the script to fit your desired SUMMIT conditions and testing files:
 
```bash
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
bash RNAseq_analyzer_mouse_180706.sh ../04_testing/metadata_mouse.txt 12
     # modify the SECOND argument to point to YOUR metadata.file
     # modify the THIRD argument to indicate the number of THREADS you 
     # want to use. This number must match the number in #SBATCH --ntasks=#

## clean up by zipping .fastq files and deleting extra files
#bash RNAseq_cleanup_mouse_180706.sh ../04_testing/metadata_mouse.txt
     # modify the SECOND argument to point to YOUR metadata.file
```
* First, double check that the inputdir specified within the program ```RNAseq_cleanup_mouse_180706.sh``` points to ```../04_testing/```.
* Double check that the analyzer will run first. Do this by ensuring the command ```bash RNAseq_analyzer_mouse_180706.sh ../04_testing/metadata_mouse.txt 12``` is uncommented (has no # in front of it). Ensure that the next command ```#bash RNAseq_cleanup_mouse_180706.sh ../04_testing/metadata_mouse.txt``` IS commented (does have a # sign in front of it).
* Run the ```execute_RNAseq_pipeline.sh``` script by entering...
```
$ sbatch execute_RNAseq_pipeline.sh
```
* Next, run the cleanup function but commenting the command ```#bash RNAseq_analyzer_mouse_180706.sh ../04_testing/metadata_mouse.txt 12``` and un-commenting the second: ```bash RNAseq_cleanup_mouse_180706.sh ../04_testing/metadata_mouse.txt```. Then execute:
```
$ sbatch execute_RNAseq_pipeline.sh
```
* Check that the resulting 03_output directory matches the example output directory in ```.../04_testing/2018-07-09_output/"


 ## 4: Modify wrapper... execute_RNAseq_pipeline.sh 
 * Open **execute_RNAseq_pipeline.sh** and modify the script to fit YOUR desired SUMMIT conditions, intput directories, input files, and other preferences.
 
```
#SBATCH --job-name=execute_RNAseq_pipeline 
#SBATCH --nodes=1
#SBATCH --ntasks=24      # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --partition=shas
#SBATCH --qos=normal     # modify this to reflect which queue you want to use. Options are 'normal' and 'testing'
#SBATCH --time=4:00:00   # modify this to reflect how long to let the job go. This indicates 4 hours.
#SBATCH --output=log_RNAseq_pipe2_%j.txt

##
source /projects/dcking@colostate.edu/paths.bashrc

## execute the RNA-seq_pipeline
bash RNAseq_analyzer_mouse_180706.sh ../01_input/metadata_mouse.txt 24  
     # modify the SECOND argument to point to YOUR metadata.file
     # modify the THIRD argument to indicate the number of THREADS you 
     # want to use. This number must match the number in #SBATCH --ntasks=#

## clean up by zipping .fastq files and deleting extra files
bash RNAseq_cleanup_mouse_180706.sh ../01_input/metadata_mouse.txt
     # modify the SECOND argument to point to YOUR metadata.file
```
* Double check that both ``` RNAseq_analyzer_mouse_180706.sh``` and ```RNAseq_cleanup_mouse_180706.sh``` meet your desired input files, input directories, file structures, and preferences.
* Execute in either in two steps (analyzer first, then clean up second, as in the testing phase, above) or as one step. 
* Execute as follows...
```
$ sbatch execute_RNAseq_pipeline.sh
```
* Save the resulting log files.
