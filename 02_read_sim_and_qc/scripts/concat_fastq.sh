#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - small queue      |
#   |    - 1 CPU + 1 Gb     |
#   +-----------------------+
#
#  Replace the USER name in this script with your username and
#  call your project whatever you want
#
#  This script must be made executable like this
#    chmod +x my_script
#
#  Submit this script to the queue with a command like this
#    run_script my_script.sh
#
#  My preferred setup before running:
#    -- script to be run in /home/projectdir/scripts
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
# LIST_FILE_NAME=sample_cat_test.txt

## Set input directory where fastq files to concatenate are stored
INPUT_DIR=/scratch/${USER}_${PROJ}/input_sim_fastq_files/

## Set output directory - where to save concatenated fastq files
OUTPUT_DIR_NAME=output_sim_fastq_files_concat
OUTPUT_DIR=/scratch/${USER}_${PROJ}/${OUTPUT_DIR_NAME}/

# Assume project directory has already been created
# /scratch/${USER}_${PROJ}/

# cd to project directory
cd /scratch/${USER}_${PROJ}/

# create directory to save concatenated fastq files
mkdir ${OUTPUT_DIR_NAME}
chmod 700 ${OUTPUT_DIR_NAME}

# ## cd into directory with output files from art
cd ${INPUT_DIR}

# Concatenate files 

while read -a line
	do
        cat ${line[0]}_11.fq ${line[0]}_21.fq > ${OUTPUT_DIR}${line[0]}_f.fq
        cat ${line[0]}_12.fq ${line[0]}_22.fq > ${OUTPUT_DIR}${line[0]}_r.fq
  	done < ${LIST_DIR}${LIST_FILE_NAME}
 mail -s 'FASTQ concatenate finished' kirkseykb1@appstate.edu <<< 'Concatenation finished'
