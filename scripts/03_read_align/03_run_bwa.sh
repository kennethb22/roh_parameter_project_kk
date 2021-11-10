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
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=03_read_align
PREV_STEP=02_read_sim_and_qc
SCRIPT=03_run_bwa.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

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

cd ${OUTPUT_DIR}

while read -a line; do

    OUT_FILE=${line[0]}_genome.bam
    start_logging "bwa align - ${OUT_FILE}"

    bwa mem -t 20 -M \
        ${REF_GENOME_FILE} \
        ${INPUT_DIR}/${line[0]}_f.fq ${INPUT_DIR}/${line[0]}_r.fq >${OUT_FILE}

    stop_logging

done <${SAMPLE_ID_LIST}

# -----------------------------------------------------------------------------
# Sort aligned read bam files
# -----------------------------------------------------------------------------

while read -a line; do

    OUT_FILE=${line[0]}_genome_sorted.bam
    start_logging "samtools sort - ${OUT_FILE}"

    samtools sort -@ 19 \
        -o ${OUT_FILE} \
        ${line[0]}_genome.bam

    stop_logging

done <${SAMPLE_ID_LIST}

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
