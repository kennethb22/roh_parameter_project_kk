#!/bin/bash
#
#   +-------------------------+
#   |  USE:                   |
#   |    - SMALL queue        |
#   |    - 5 CPU + def Gb     |
#   +-------------------------+
#
#  Replace the USER name in this script with your username and
#  call your project whatever you want
#
#  This script must be made executable like this
#    chmod +x my_script
#
#  Submit this script to the queue with a command like this
#    run_script_scratch my_script.sh
#
#  My preferred setup before running:
#    -- script to be run in /home/scripts
#    -- project directory (of same name as script) in /home/
#    -- /input/ and /output/ subdirs within project dir

##  Set username
USER=aubkbk001

## Set project name
PROJ=02_read_sim_and_qc

## Set directory where list of sample files is stored
LIST_DIR=/scratch/aubkbk001_01_slim/

## Set file name of file where list of sample files is stored
LIST_FILE_NAME=sample_id_list_m5e-07_r1e-8_p500.txt

## Set input directory where fastq files to concatenate are stored
INPUT_DIR=/scratch/${USER}_${PROJ}/output_sim_fastq_files_concat/

## Set output directory - where to save concatenated fastq files
OUTPUT_DIR_NAME=output_rawread_fastqc
OUTPUT_DIR=/scratch/${USER}_${PROJ}/${OUTPUT_DIR_NAME}/

# Assume project directory has already been created
# /scratch/${USER}_${PROJ}/

# cd to project directory
cd /scratch/${USER}_${PROJ}/

# create directory to save concatenated fastq files
mkdir ${OUTPUT_DIR_NAME}
chmod 700 ${OUTPUT_DIR_NAME}

# ## cd into directory with concatenated fastq files from art
cd ${INPUT_DIR}

## --------------------------------
## Load modules 
# module load trimmomatic/0.38
module load fastqc/0.11.9

## Run fastqc on 2 read files per sample

while read -a line
do
	fastqc -t 5 -o ${OUTPUT_DIR} \
	${line[0]}_f.fq \
	${line[0]}_r.fq 

done < ${LIST_DIR}${LIST_FILE_NAME}

 mail -s 'FASTQC run finished' kirkseykb1@appstate.edu <<< 'FASTQC run finished'
