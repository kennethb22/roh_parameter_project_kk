#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - LARGE queue      |
#   |    - 20 CPU + 20 Gb   |
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
#    -- script to be run in /home/scripts
#    -- project directory (of same name as script) in /home/
#    -- /input/ and /output/ subdirs within project dir

# -----------------------------------------------------------------------------
# Set Variables
# -----------------------------------------------------------------------------

## Load variables from settings file
source /home/aubkbk001/roh_param_project/init_script_vars.sh

## Set location of file contaning list of sample files to process
SAMPLE_LIST_DIR=/scratch/aubkbk001_01_slim/
SAMPLE_LIST_FILE_NAME=sample_id_list_m5e-07_r1e-8_p500.txt
# SAMPLE_LIST_FILE_NAME=sample_cat_test.txt

## Set location of output directory for unsorted bam files
OUTPUT_DIR_UNSORTED=/scratch/${USER}_${PROJ}/output/unsorted_bam/

## Set location of output directory for unsorted bam files
OUTPUT_DIR_SORTED=/scratch/${USER}_${PROJ}/output/sorted_bam/

# -----------------------------------------------------------------------------
# Create working directories
# -----------------------------------------------------------------------------

## Create a directory on /scratch
mkdir /scratch/${USER}_${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/

## create output directory and subdirectories

mkdir output
chmod 700 output
mkdir ./output/sorted_bam
chmod 700 ./output/sorted_bam
mkdir ./output/unsorted_bam
chmod 700 ./output/unsorted_bam

## create input directory
mkdir input
chmod 700 input

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load bwa/0.7.12
module load samtools/1.11

# -----------------------------------------------------------------------------
# Index reference genomes
# -----------------------------------------------------------------------------

bwa index ${REF_GENOME_FILE}

# -----------------------------------------------------------------------------
# Align reads to reference genome
# -----------------------------------------------------------------------------

while read -a line; do
    bwa mem -t 20 -M \
        ${REF_GENOME_FILE} \
        ./input/${line[0]}_f.fq ./input/${line[0]}_r.fq >${OUTPUT_DIR_UNSORTED}${line[0]}_genome.bam
done <${SAMPLE_LIST_DIR}${SAMPLE_LIST_FILE_NAME}

# -----------------------------------------------------------------------------
# Sort aligned read bam files
# -----------------------------------------------------------------------------

while read -a line; do
    samtools sort -@ 19 \
        -o ${OUTPUT_DIR_SORTED}${line[0]}_genome_sorted.bam \
        ${OUTPUT_DIR_UNSORTED}${line[0]}_genome.bam
done <${SAMPLE_LIST_DIR}${SAMPLE_LIST_FILE_NAME}
