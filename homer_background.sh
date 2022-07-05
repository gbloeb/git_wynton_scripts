#!/bin/env bash
#!/bin/bash                         #-- what is the language of this shell

#                                  #-- Any line that starts with #$ is an instruction to SGE

#$ -S /bin/bash                     #-- the shell for the job

#$ -o ~/log                        #-- output directory (fill in)

#$ -e ~/log                        #-- error directory (fill in)

#$ -cwd                            #-- tell the job that it should start in your working directory

#$ -r y                            #-- tell the system that if a job crashes, it should be restarted

#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined

#$ -l mem_free=10G

#$ -l scratch=20G

#$ -l h_rt=6:00:00

#$ -pe smp 4

#$ -m ea                           #--email when done

#$ -M gabriel.loeb@ucsf.edu        #--email

OUTPUT_DIRECTORY=$1
INPUTBED=$2
INPUTBEDFILE="${INPUTBED##*/}"
GENOME=$3
SIZE=$4
BCKGRND=$5
NUMMOTIFS=$6


findMotifsGenome.pl $INPUTBED $GENOME "$OUTPUT_DIRECTORY"/"${INPUTBEDFILE%.*}"_"$SIZE"bp \
-size $SIZE -p "${NSLOTS:-1}" -bg $BCKGRND  -S $NUMMOTIFS

