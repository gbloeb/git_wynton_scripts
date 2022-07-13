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

#$ -l h_rt=2:00:00

#$ -pe smp 1

#$ -m ea                           #--email when done

#$ -M gabriel.loeb@ucsf.edu        #--email

OUTPUT_DIRECTORY=$1
INPUT_DIRECTORY=$2
PEAKS=$3

PEAKS_BASE="${PEAKS##*/}"

cd $OUTPUT_DIRECTORY
for FULL_SAMPLE in "$INPUT_DIRECTORY"/*_S*
do
	SAMPLE=$(basename "$FULL_SAMPLE" .bed)
	echo $SAMPLE
	READS="$INPUT_DIRECTORY"/"$SAMPLE"/aligned_mm10_exact/"$SAMPLE"_shift.bed

  	awk '$1 !~ /_/' $READS > $OUTPUT_DIRECTORY/int_files/"$SAMPLE"_stdChr.bed #Remove nonstandard chromosomes
	bedtools coverage -a $PEAKS -b $OUTPUT_DIRECTORY/int_files/"$SAMPLE"_stdChr.bed > $OUTPUT_DIRECTORY/"$SAMPLE"_coverage_"$PEAKS_BASE"
done
