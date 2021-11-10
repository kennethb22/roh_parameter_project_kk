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

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=02_read_sim_and_qc
PREV_STEP=01_run_slim
SCRIPT=02b_concat_fastq.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Concatenate forward and reverse read files to a single file for each
# individual
#
# OUTPUT_DIR and SAMPLE_ID_LIST are set in init_script_vars.sh
# -----------------------------------------------------------------------------

# cd to project directory
cd ${OUTPUT_DIR}

# Concatenate files

while read -a line; do

    start_logging "Concat fastq - ${line[0]}"

    cat ${line[0]}_11.fq ${line[0]}_21.fq >${line[0]}_f.fq
    cat ${line[0]}_12.fq ${line[0]}_22.fq >${line[0]}_r.fq

    stop_logging

done <${SAMPLE_ID_LIST}

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
