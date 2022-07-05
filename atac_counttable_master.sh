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
OUTPUT_DIRECTORY=$2 #Will make a directory within this with the correct qval and lambda of called peaks
QVAL=$3
NOLAMBDA=$4

module load CBI bedtools2 r

mkdir $OUTPUT_DIRECTORY

#Now reset OUTPUT_DIRECTORY to correspond to the subdirectory

if [ "$NOLAMBDA" = true ]
then
OUTPUT_DIRECTORY=$OUTPUT_DIRECTORY/nolambda_q"$QVAL"
PEAKS=comb_shift.bed_nolambda_q"$QVAL"_noshift_peaks.narrowPeak
else
OUTPUT_DIRECTORY=$OUTPUT_DIRECTORY/wlambda_q"$QVAL"
PEAKS=comb_shift.bed_wlambda_q"$QVAL"_noshift_peaks.narrowPeak
fi

PEAKS_BASE="${PEAKS##*/}"





mkdir $OUTPUT_DIRECTORY
mkdir $OUTPUT_DIRECTORY/int_files


#Make bed file from all the samples to call peaks on
cat $INPUT_DIRECTORY/*/aligned_mm10_exact/*_shift.bed > $OUTPUT_DIRECTORY/int_files/comb_shift.bed
ls $INPUT_DIRECTORY/*/aligned_mm10_exact/*_shift.bed > $OUTPUT_DIRECTORY/int_files/files_in_comb_shift_bed.txt

#Call peaks
job_id=$(qsub -terse ~/git_wynton_scripts/macs3_callpeaks_frombed.sh \
$OUTPUT_DIRECTORY/int_files/comb_shift.bed \
$OUTPUT_DIRECTORY/int_files/ \
2.3e9 \
$QVAL \
$NOLAMBDA)




#Calculate peak coverage


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
Rscript ~/git_wynton_scripts/220703_gen_countTable_calldiffpeaks.R  $OUTPUT_DIRECTORY

mkdir $OUTPUT_DIRECTORY/homer_output

#Run homer
#given bp
qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.05.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.01.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.001.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.05.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

#500bp
qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.05.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.01.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.001.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.05.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

#250bp
qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.05.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.01.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.001.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub ~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.05.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10
