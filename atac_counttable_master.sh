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

#$ -l h_rt=00:29:30

#$ -pe smp 1

#$ -m ea                           #--email when done

#$ -M gabriel.loeb@ucsf.edu        #--email



INPUT_DIRECTORY=$1
OUTPUT_DIRECTORY=$2 #Will make a directory within this with the correct qval and lambda of called peaks
QVAL=$3
NOLAMBDA=$4
CTRL1=$5
CTRL2=$6
CTRL3=$7
CTRL4=$8
CTRL5=$9

module load CBI bedtools2 r

echo "CTRL1"
echo $CTRL1

echo "CTRL2"
echo $CTRL2

echo "CTRL3"
echo $CTRL3

echo "CTRL4"
echo $CTRL4

echo "CTRL5"
echo $CTRL5



mkdir $OUTPUT_DIRECTORY

#Now reset OUTPUT_DIRECTORY to correspond to the subdirectory

if [ "$NOLAMBDA" = true ]
then
OUTPUT_DIRECTORY=$OUTPUT_DIRECTORY/nolambda_q"$QVAL"
PEAKS=$OUTPUT_DIRECTORY/int_files/comb_shift.bed_nolambda_q"$QVAL"_noshift_peaks.narrowPeak
else
OUTPUT_DIRECTORY=$OUTPUT_DIRECTORY/wlambda_q"$QVAL"
PEAKS=$OUTPUT_DIRECTORY/int_files/comb_shift.bed_wlambda_q"$QVAL"_noshift_peaks.narrowPeak
fi



mkdir $OUTPUT_DIRECTORY
mkdir $OUTPUT_DIRECTORY/int_files
mkdir $OUTPUT_DIRECTORY/log_files

#Make bed file from all the samples to call peaks on
cat $INPUT_DIRECTORY/*/aligned_mm10_exact/*_shift.bed > $OUTPUT_DIRECTORY/int_files/comb_shift.bed
ls $INPUT_DIRECTORY/*/aligned_mm10_exact/*_shift.bed > $OUTPUT_DIRECTORY/int_files/files_in_comb_shift_bed.txt

#Call peaks
qsub -N callpeaks_"$QVAL"_"$NOLAMBDA" -o $OUTPUT_DIRECTORY/log_files/callpeaks_"$QVAL"_"$NOLAMBDA".log -j y  \
~/git_wynton_scripts/macs3_callpeaks_frombed_noshift.sh \
$OUTPUT_DIRECTORY/int_files/comb_shift.bed \
$OUTPUT_DIRECTORY/int_files/ \
2.3e9 \
$QVAL \
$NOLAMBDA



#Calculate peak coverage

qsub -N peakcov_"$QVAL"_"$NOLAMBDA" -hold_jid callpeaks_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/peakcov_"$QVAL"_"$NOLAMBDA".log -j y  \
~/git_wynton_scripts/atac_counttable_bedcoverage.sh \
$OUTPUT_DIRECTORY \
$INPUT_DIRECTORY \
$PEAKS


#Make count table, homer input files, volcano plot
qsub -N count_table_"$QVAL"_"$NOLAMBDA" -hold_jid peakcov_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/count_table_"$QVAL"_"$NOLAMBDA".log -j y  \
~/git_wynton_scripts/220703_gen_countTable_calldiffpeaks.sh  $OUTPUT_DIRECTORY $CTRL1 $CTRL2 $CTRL3 $CTRL4 $CTRL5 

mkdir $OUTPUT_DIRECTORY/homer_output

#Run homer
#given bp
qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.05_given.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.05.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.01_given.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.01.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.001_given.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.001.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.05_given.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.05.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.01_given.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.01.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.001_given.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.001.bed \
mm10 \
given \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

#500bp
qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.05_500.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.05.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.01_500.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.01.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.001_500.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.001.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.05_500.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.05.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.01_500.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.01.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.001_500.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.001.bed \
mm10 \
500 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

#250bp
qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.05_250.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.05.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.01_250.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.01.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_up_0.001_250.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_up_0.001.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.05_250.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.05.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.01_250.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.01.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10

qsub -hold_jid count_table_"$QVAL"_"$NOLAMBDA" \
-o $OUTPUT_DIRECTORY/log_files/homer_background_down_0.001_250.log -j y  \
~/git_wynton_scripts/homer_background.sh $OUTPUT_DIRECTORY/homer_output \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_down_0.001.bed \
mm10 \
250 \
$OUTPUT_DIRECTORY/homer_input/TestedPeaks_background.bed \
10