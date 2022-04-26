#!/bin/bash
#SBATCH --job-name=02_read_sim_and_qc
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 02:00:00
#SBATCH --mem=4000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu

#  This script must be made executable like this
#    chmod +x my_script

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=02_read_sim_and_qc
PREV_STEP=01_run_slim
SCRIPT=02c_fastqc.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# Create directory to save fastqc reports
cd ${OUTPUT_DIR}
cd ..
FASTQC_OUT_DIR=fastqc_reports
mkdir ${FASTQC_OUT_DIR}

## --------------------------------
## Load modules# module load trimmomatic/0.38
module load fastqc/0.11.9

## Run fastqc on 2 read files per sample

while read -a line; do

    start_logging "fastqc - ${line[0]}"

    fastqc -t 5 -o ./${FASTQC_OUT_DIR} \
        ${OUTPUT_DIR}/${line[0]}_f.fq \
        ${OUTPUT_DIR}/${line[0]}_r.fq

    stop_logging

done <${SAMPLE_ID_LIST}

# mail -s 'FASTQC run finished' kirkseykb1@appstate.edu <<<'FASTQC run finished'
