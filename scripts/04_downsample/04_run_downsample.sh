#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - ?LARGE queue     |
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

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=04_downsample
PREV_STEP=03_read_align
SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load samtools/1.11

# -----------------------------------------------------------------------------
# Do the downsampling.
# -----------------------------------------------------------------------------

for i in $(seq 0 $cvgCnt); do

    # Create output directory for each coverage level.

    CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_cvg_${cvgX[i]}
    mkdir ${CVG_OUTPUT_DIR}

    # Downsample each individual

    while read -a line; do

        OUT_FILE=${line[0]}_cvg_${cvgX[i]}.bam
        start_logging "samtools view - ${OUT_FILE}"

        samtools view -s ${cvgP[i]} -@ 19 \
            -o ${CVG_OUTPUT_DIR}/${OUT_FILE} \
            ${INPUT_DIR}/${line[0]}_genome_sorted.bam

        stop_logging

    done <${SAMPLE_ID_LIST}
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
