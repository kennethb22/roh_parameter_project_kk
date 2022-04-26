#!/bin/bash
#
#SBATCH --job-name=02b_concat_fastq
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mem=4000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=02_read_sim_and_qc
PREV_STEP=01_run_slim
SCRIPT=02b_concat_fastq.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

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

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
