#!/bin/env bash
#!/bin/bash                         #-- what is the language of this shell

#                                  #-- Any line that starts with #$ is an instruction to SGE

#$ -S /bin/bash                     #-- the shell for the job

#$ -o ~/log                        #-- output directory (fill in)

#$ -e ~/log                        #-- error directory (fill in)

#$ -cwd                            #-- tell the job that it should start in your working directory

#$ -r y                            #-- tell the system that if a job crashes, it should be restarted

#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined

#$ -l mem_free=50G

#$ -l scratch=50G

#$ -l h_rt=20:29:30

#$ -pe smp 1

#$ -m ea                           #--email when done

#$ -M gabriel.loeb@ucsf.edu        #--email



INPUT_DIRECTORY=$1
PEAKS=$2
OUTPUT_DIRECTORY=$3

PEAKS_BASE="${PEAKS##*/}"

module load CBI bedtools2 r


#Call peaks

#Calculate peak coverage
mkdir $OUTPUT_DIRECTORY
mkdir $OUTPUT_DIRECTORY/int_files

cd $OUTPUT_DIRECTORY

for FULL_SAMPLE in "$INPUT_DIRECTORY"/*_S*
do
	SAMPLE=$(basename "$FULL_SAMPLE" .bed)
	echo $SAMPLE
	READS="$INPUT_DIRECTORY"/"$SAMPLE"/aligned_mm10_exact/"$SAMPLE"_shift.bed

  	awk '$1 !~ /_/' $READS > $OUTPUT_DIRECTORY/int_files/"$SAMPLE"_stdChr.bed #Remove nonstandard chromosomes
	bedtools coverage -a $PEAKS -b $OUTPUT_DIRECTORY/int_files/"$SAMPLE"_stdChr.bed > $OUTPUT_DIRECTORY/"$SAMPLE"_coverage_"$PEAKS_BASE"
done

#Make count table, homer input files, volcano plot
Rscript 220703_gen_countTable_calldiffpeaks.R $OUTPUT_DIRECTORY
