#!/bin/bash
#
#SBATCH --job-name=03_downsample
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 04:00:00
#SBATCH --mem=20000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu

# srun -N1 -n20 -t02:00:00 --mem=20000 --pty bash

#   +-----------------------+
#   |  USE:                 |
#   |    - ?LARGE queue     |
#   |    - 20 CPU + 20 Gb   |
#   +-----------------------+
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=04_downsample
PREV_STEP=03_read_align
# This doesn't work when running via sbatch on easley. It does work when
# running the script in the shell manually, or running it under srun. Go figure.
# SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")
SCRIPT=04_run_downsample.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load samtools/1.11

# -----------------------------------------------------------------------------
# Do the downsampling.
# -----------------------------------------------------------------------------

for i in $(seq 0 $cvgCnt); do

    # # Create output directory for each coverage level.

    CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_cvg_${cvgX[i]}
    mkdir ${CVG_OUTPUT_DIR}

    # # Downsample each individual

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

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
