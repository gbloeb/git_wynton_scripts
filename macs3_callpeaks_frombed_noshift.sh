
#!/bin/env bash
#!/bin/bash                         #-- what is the language of this shell

#                                  #-- Any line that starts with #$ is an instruction to SGE

#$ -S /bin/bash                     #-- the shell for the job

#$ -o ~/log                        #-- output directory (fill in)

#$ -e ~/log                        #-- error directory (fill in)

#$ -cwd                            #-- tell the job that it should start in your working directory

#$ -r y                            #-- tell the system that if a job crashes, it should be restarted

#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined

#$ -l mem_free=20G

#$ -l scratch=200G

#$ -l h_rt=10:00:00

#$ -pe smp 1

#$ -m ea                           #--email when done

#$ -M gabriel.loeb@ucsf.edu        #--email

INPUT=$1
OUTPUTDIRECTORY=$2
GS=$3 # 2.3e9 for mm10 with 50bp PE reads
QVAL=$4
NOLAMBDA=$5


INPUTFILENAME="${INPUT##*/}"


conda activate macs3

if [ "$NOLAMBDA" = true ]
then

macs3 callpeak -t $INPUT \
-f BED --outdir $OUTPUTDIRECTORY -n "$INPUTFILENAME"_nolambda_q"$QVAL"_noshift \
 -g $GS --nomodel --keep-dup all -q $QVAL --nolambda

#callpeaks with lambda
else

macs3 callpeak -t $INPUT \
-f BED --outdir $OUTPUTDIRECTORY -n "$INPUTFILENAME"_wlambda_q"$QVAL"_noshift \
 -g $GS  --nomodel --keep-dup all -q $QVAL

fi




