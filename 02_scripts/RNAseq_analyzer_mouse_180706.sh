#!/usr/bin/bash/

################################################
# Program:
# RNAseq_analyzer_mouse.sh
#
# Author:
# Erin Osborne Nishimura
#
# Date initiated:
# July 3, 2018
#
# Description:
#
# This is a very basic RNA-seq pipeline that I use for analyzing mouse RNA-seq paired-end sequencing
#
# Requires: fastqc, hisat2, htseq, samtools
#
# Requires: For each sample, forward and reverse paired end sequences are required. The .fastq.gz files
#           should be placed in the inputdir directory.
#
# Requires: A metadata file with two columns. The first two columns are fastq.gz file names. 
#           The third column is a "nickname" of each sample
#			Later columns can be included with other metadata information
#           Metadata file should be placed within the inputdir directory.
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



echo -e ">>> INITIATING analyzer with command:\n\t$0 $@"


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





#This is the output_directory:
DATE=`date +%Y-%m-%d`
outputdir="../03_output/"$DATE"_output/"


echo -e ">>> MAKING output directory"
echo -e "\tmkdir $outputdir"
mkdir -p $outputdir



####### META DATA #############



#These are the sample names, R1:
samples1=( $(cut -f 1 --output-delimiter=' ' $1) )

#These are the sample names, R2:
samples2=( $(cut -f 2 --output-delimiter=' ' $1) )

#These are the nicknames I want to give the files:
names=( $(cut -f 3 --output-delimiter=' ' $1) )



####### PIPELINE ##############

# Report back to the user which files will be processed and which names they'll be given:
echo -e ">>> INPUT: This script will process files from the metafile:\n\t$1"
echo -e ">>> PLAN: This script will process the sample files into the following names: "
echo -e "\tSAMPLE1\tSAMPLE2\tNAMES"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    echo -e "\t${samples1[$counter]}\t${samples2[$counter]}\t${names[$counter]}"
done



# UNZIP all the input files: 
echo -e "\n>>> GUNZIP: unzipping _R1.fastq.gz files:"
for fastqfile in ${samples1[@]}
do
   echo -e "\tgunzipping $fastqfile"
   #gunzip $inputdir$fastqfile
done

echo -e ">>> GUNZIP: unzipping _R2.fastq.gz files:"
for fastqfile in ${samples2[@]}
do
   echo -e "\tgunzipping $fastqfile"
   #gunzip $inputdir$fastqfile
done



# FASTQC to determine quality
echo -e "\n>>> FASTQC: analyzing quality of each .fastq file"
mkdir -p $outputdir"01_fastqc"

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	samplename=${names[$counter]}
	sample1=${samples1[$counter]}
	sample2=${samples2[$counter]}
	
	#Chop off the .gz of each name:
	unzippedfile1=${sample1//.gz/}
	unzippedfile2=${sample2//.gz/}

	#Make output directories
   	mkdir -p $outputdir"01_fastqc/"$samplename"_1"
   	mkdir -p $outputdir"01_fastqc/"$samplename"_2"	
   	
   	# execute fastqc
   	echo -e "\t$ fastqc -o $outputdir"01_fastqc/"$samplename"_1" -t 20 $inputdir$unzippedfile1"
   	#time fastqc -o $outputdir"01_fastqc/"$samplename"_1" -t 20 $inputdir$unzippedfile1
   	
   	echo -e "\t$ fastqc -o $outputdir"01_fastqc/"$samplename"_2" -t 20 $inputdir$unzippedfile2"
   	#time fastqc -o $outputdir"01_fastqc/"$samplename"_2" -t 20 $inputdir$unzippedfile2
	
done



# HISAT2 to align to the genome
echo -e "\n>>> HISAT2: aligning each sample to the genome"
outhisat2=$outputdir"02_hisat2/"
mkdir -p $outhisat2

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
	samplename=${names[$counter]}
	sample1=${samples1[$counter]}
	sample2=${samples2[$counter]}
	
	#Chop off the .gz of each name:
	unzippedfile1=${sample1//.gz/}
	unzippedfile2=${sample2//.gz/}

   	# execute hisat2
   	echo -e "\t$ hisat2 -x $hisat2path -1 ${inputdir}${unzippedfile1} -2 ${inputdir}${unzippedfile2} -S ${outhisat2}${samplename}.sam --summary-file ${outhisat2}${samplename}_summary.txt --no-unal -p $pthread"
	#time hisat2 -x $hisat2path -1 ${inputdir}${unzippedfile1} -2 ${inputdir}${unzippedfile2} -S ${outhisat2}${samplename}.sam --summary-file ${outhisat2}${samplename}_summary.txt --no-unal -p $pthread

done



# RUN HT-SEQ ON ALL FILES USING ORGANISM GTF FILE:
echo -e "\n>>> HTSEQ: Run HTSeq-counts on all files to tabulate read counts per gene"
outhtseq=$outputdir"04_htseq/"
mkdir -p $outhtseq

for (( counter=0; counter < ${#samples1[@]}; counter++ ))
do
    samfile=${names[$counter]}.sam
    echo -e $samfile
    #samtools view --threads 19 -h -o ${outhtseq}${samfile} ${samout}${bamfile}
    
    echo -e "\thtseq-count --stranded=yes --minaqual=20 -r pos -q ${outhtseq}${samfile} $gtffile > ${outhtseq}${names[$counter]}_counts.txt"
    #time htseq-count --stranded=yes --minaqual=20 -r pos -q ${outhisat2}${samfile} $gtffile > ${outhtseq}${names[$counter]}_counts.txt


    #rm ${outhtseq}${samfile}
done



# SAMTOOLS and BAMCOVERAGE: to convert .sam output to uploadable .bam and .wg files
echo -e "\n>>> SAMTOOLS/BAMCOVERAGE: to convert files to uploadable _sort.bam and _sort.bam.bai files:"
samout=$outputdir"03_samtools/"
mkdir -p $samout

for seqname in ${names[@]}
do
    # echo
    echo -e "Samtools convert: ${seqname}"
    
    # Samtools: compress .sam -> .bam
	echo -e "\t$ samtools view --threads $pthread -bS ${outhisat2}${seqname}.sam > ${samout}${seqname}.bam"
	samtools view --threads $pthread -bS ${outhisat2}${seqname}.sam > ${samout}${seqname}.bam
    
    # Samtools: sort .bam -> _sort.bam
    echo -e "\t$ samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam"
    #samtools sort --threads $pthread -o ${samout}${seqname}_sort.bam --reference $genomefa ${samout}${seqname}.bam
    
    # Samtools: index _sort.bam -> _sort.bam.bai
    echo -e "\t$ samtools index ${samout}${seqname}_sort.bam"
    #samtools index ${samout}${seqname}_sort.bam
    
    # bamCoverage: 
    export PYTHONPATH=/projects/dcking@colostate.edu/lib/python3.5/site-packages/
    echo -e "BamCoverage convert: ${seqname}"
    echo -e "\t$ bamCoverage -b ${samout}${seqname}_sort.bam -o ${samout}${seqname}_sort.bw --outFileFormat bigwig -p $pthread --normalizeUsing CPM --binSize 1"
    time bamCoverage -b ${samout}${seqname}_sort.bam -o ${samout}${seqname}_sort.bw --outFileFormat bigwig -p $pthread --normalizeUsing CPM --binSize 1
    export PYTHONPATH=/projects/dcking@colostate.edu/lib/python2.7/site-packages/
done




##############################################
#                ERCC                        #
##############################################

echo -e "ERCC variable is $ercc"

if [ $ercc == "TRUE" ]
then
	# HISAT2 to align to the ERCC spike-in controls
	echo -e "\n>>> HISAT2: aligning each sample to the ERCC spike in library"
	outercc=$outputdir"05_ercc/"
	mkdir -p $outercc

	for (( counter=0; counter < ${#samples1[@]}; counter++ ))
	do
		samplename=${names[$counter]}
		sample1=${samples1[$counter]}
		sample2=${samples2[$counter]}
	
		#Chop off the .gz of each name:
		unzippedfile1=${sample1//.gz/}
		unzippedfile2=${sample2//.gz/}

		# execute hisat2
		echo -e "\t$ hisat2 -x $erccpath -1 $unzippedfile1 -2 $unzippedfile2 -S ${outercc}${samplename}_ercc.sam --summary-file ${outercc}${samplename}_ercc_summary.txt --no-unal -p 20"
		#time hisat2 -x $erccpath -1 ${inputdir}${unzippedfile1} -2 ${inputdir}${unzippedfile2} -S ${outercc}${samplename}_ercc.sam --summary-file ${outercc}${samplename}_ercc_summary.txt --no-unal -p 20

	done


	# RUN HT-SEQ ON ALL FILES USING ERCC GTF FILE:
	echo -e "\n>>> HTSEQ_ERCC: Run HTSeq-counts on ERCC 'Genes'"
	out_ercc_counts=$outputdir"06_ercc_counts/"
	mkdir -p $out_ercc_counts

	for (( counter=0; counter < ${#samples1[@]}; counter++ ))
	do

		echo -e "samtools sort --output-fmt SAM -o ${out_ercc_counts}${names[$counter]}_ercc_sort.sam -n ${outercc}${names[$counter]}_ercc.sam"
		samtools sort --output-fmt SAM -o ${out_ercc_counts}${names[$counter]}_ercc_sort.sam -n ${outercc}${names[$counter]}_ercc.sam
	
		samfile=${names[$counter]}_ercc_sort.sam
		echo -e $samfile
	
		echo -e "\thtseq-count --stranded=yes --minaqual=20 -r pos -q ${out_ercc_counts}${samfile} $erccgtf > ${out_ercc_counts}${names[$counter]}_counts.txt"
		time htseq-count --stranded=yes --minaqual=20 -r pos -q ${out_ercc_counts}${samfile} $erccgtf > ${out_ercc_counts}${names[$counter]}_counts.txt

		rm ${out_ercc_counts}${samfile}
	done
	
fi



######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> FASEQC VERSION:"
fastqc --version
echo -e "\n>>> HISAT2 VERSION:"
hisat2 --version
echo -e "\n>>> SAMTOOLS VERSION:"
samtools --version
echo -e "\n>>> HTSEQ-COUNT VERSION:"
htseq-count --help | grep "version"
echo -e "\n>>> BAMCOVERAGE VERSION:"
export PYTHONPATH=/projects/dcking@colostate.edu/lib/python3.5/site-packages/
bamCoverage --version
echo -e ">>> END: Analayzer complete."
