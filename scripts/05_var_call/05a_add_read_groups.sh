#!/bin/bash
#
#  +---------+
#  | BATCHED |
#  +---------+
#
#  +----------------------+
#  | REQUEST 4 CPU + 60gb |
#  +----------------------+
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
#  In user home directory, make a backup of the .asc_queue file and replace it
#  with .asc_queue_gatk_hc file before running haplotype calling. The .asc_queue
#  file sets the default parameters used by run_script.
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
sleep 5

module load picard/1.79
# module load gatk/4.1.4.0

# -----------------------------------------------------------------------------
# Run Haplotype caller on all sample files
# -----------------------------------------------------------------------------

# Iterate over levels of coverage ------------------------------------------
declare -a cvgX=(30x)
cvgCnt=${#cvgX[@]}
let cvgCnt-=1

for i in $(seq 0 $cvgCnt); do

    # Create output directory for each coverage level.

    CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_cvg_${cvgX[i]}
    mkdir ${CVG_OUTPUT_DIR}

    # Set input directory for this coverage

    CVG_INPUT_DIR=${INPUT_DIR}/sample_cvg_${cvgX[i]}

    # Process each individual sample file ----------------------------------

    while read -a line; do

        # # ------------------------------------------------------------------
        # # Add read group information
        # # ------------------------------------------------------------------

        IN_FILE=${line[0]}_cvg_${cvgX[i]}.bam
        OUT_FILE=${line[0]}_cvg_${cvgX[i]}_rgroups.bam
        RGROUPS_FILE=${OUT_FILE}
        start_logging "picard/AddOrReplaceReadGroups - ${OUT_FILE}"

        ## Do picard read groups

        java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar \
            I=${CVG_INPUT_DIR}/${IN_FILE} \
            O=${CVG_OUTPUT_DIR}/${OUT_FILE} SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=${line[0]} RGPU=unit4

        stop_logging

        # ------------------------------------------------------------------
        # Index BAM files
        # ------------------------------------------------------------------

        OUT_FILE=${line[0]}_cvg_${cvgX[i]}_rgroups.bai
        start_logging "picard/BuildBamIndex - ${OUT_FILE}"

        ## Run picard BuildBamIndex

        java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/BuildBamIndex.jar \
            I=${CVG_OUTPUT_DIR}/${RGROUPS_FILE}

        stop_logging

    done <${SAMPLE_ID_LIST}
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
