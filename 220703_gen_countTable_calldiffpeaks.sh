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


OUTPUT_DIRECTORY=$1
CTRL1=$2
CTRL2=$3
CTRL3=$4
CTRL4=$5
CTRL5=$6
CTRL6=$7



module load CBI r

Rscript ~/git_wynton_scripts/220703_gen_countTable_calldiffpeaks.R  $OUTPUT_DIRECTORY $CTRL1 $CTRL2 $CTRL3 $CTRL4 $CTRL5 $CTRL6
