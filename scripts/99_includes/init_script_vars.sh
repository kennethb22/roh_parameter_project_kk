#!/bin/bash
#
# initialize_script_variables.sh
#
# This script initializes variables that are used in multiple scripts in the
# ROH calling workflow. It should be included in all scripts for each step of
# the workflow. Put the following line in the scripts:
#
#     source /home/aubkbk001/roh_param_project/init_script_vars.sh
#
# replacing aubkbk001 with your username.

##  Set username
USER=aubkbk001

## Set project name
PROJECT=roh_param_project

## Set user home directory
USER_HOME_DIR=/home/${USER}/

## Set email address
EMAIL=kirkseykb1@appstate.edu

## Set initial step name
INIT_STEP=01_slim

# -----------------------------------------------------------------------------
# Create working directories
# -----------------------------------------------------------------------------

## Create variable pointing to user's home directory for this project and step

HOME_STEP_DIR=/home/${USER}/${PROJECT}/data/${STEP}/

## Set the input directory, which is the output directory from the previous
## step
INPUT_DIR=/scratch/${USER}/${PROJECT}/data/${PREV_STEP}/output

## Create a working directory on /scratch
mkdir -p /scratch/${USER}/${PROJECT}/data/${STEP}
chmod -R 700 /scratch/${USER}
WORK_DIR=/scratch/${USER}/${PROJECT}/data/${STEP}

## Create output directory in working directory
mkdir ${WORK_DIR}/output
chmod -R 700 ${WORK_DIR}
OUTPUT_DIR=${WORK_DIR}/output

# -----------------------------------------------------------------------------
## Set location of reference genome file
# -----------------------------------------------------------------------------

REF_GENOME_FILE_PATH=/scratch/${USER}/${PROJECT}/data/${INIT_STEP}/output/slim_output_files_m5e-07_r1e-8_p500

REF_GENOME_FILE_NAME=ancestral

REF_GENOME_FILE=${REF_GENOME_FILE_PATH}/${REF_GENOME_FILE_NAME}.fasta

# -----------------------------------------------------------------------------
# Create arrays defining coverage levels and population sizes. These arrays are
# used in:
#
#    04_downsample.sh
#    05a_snp_calling.sh
#    05b_genomicsdb.sh
#    06a_run_bcftoolsROH.sh
#    06b_run_plink.sh
# -----------------------------------------------------------------------------

## Create arrays of downsample levels to be used.
## Coverage level in NNx for display in output file names
# declare -a cvgX=(50x 30x 15x 10x 05x)
# declare -a cvgX=(50x 30x)
declare -a cvgX=(05x)

## Coverage level fraction to supply to samtools
# declare -a cvgP=(1.0 0.6 0.3 0.2 0.1)
# declare -a cvgP=(1.0 0.6)
declare -a cvgP=(0.1)

## Get length of the coverage level arrays. Subtract 1 because arrays are zero
## based, and we'll iterate over the arrays from 0 to cvgCnt
cvgCnt=${#cvgX[@]}
let cvgCnt-=1

## Create array of population sizes we want to test
# declare -a popN=(100 50 30)
declare -a popN=(03)
# declare -a popN=(01)

# -----------------------------------------------------------------------------
# PLINK parameters - define arrays for each of the PLINK parameters we want to
# vary and test. Used in 06b_run_plink.sh
# -----------------------------------------------------------------------------

declare -a phwh=(0 1 2 3)   # Values for -homozyg-window-het
declare -a phwm=(5 10 50)   # Values for -homozyg-window-missing
declare -a phws=(50)        # Values for -homozyg-window-snp
declare -a phzd=(50)        # Values for -homozyg-density
declare -a phzg=(1000)      # Values for -homozyg-gap
declare -a phwt=(0.01 0.05) # Values for -homozyg-window-threshold
declare -a phzs=(25 50 100) # Values for -homozyg-snp
declare -a phzk=(10)        # Values for -homozyg-kb

# -----------------------------------------------------------------------------
# Create log file
# -----------------------------------------------------------------------------

# Set name of log file to track execution times

LOG_FILE=${OUTPUT_DIR}/${SCRIPT}_log.txt

# Delete log file if it exists

if [ -f "$LOG_FILE" ]; then
    rm ${LOG_FILE}
fi

# Write header to log file

printf "%-80s   %8s   %8s   %8s\n" "Action - Output" "Start" "End" "Duration" >${LOG_FILE}

# -----------------------------------------------------------------------------
# Define Functions
# -----------------------------------------------------------------------------

# Log file functions - used in all scripts

start_logging() {

    ACTION=$1
    START_TIME_HHMM=$(date +" %T")
    START_TIME=$(date +%s)
}

stop_logging() {

    END_TIME_HHMM=$(date +"%T")
    END_TIME=$(date +%s)
    ELAPSED=$(expr $END_TIME - $START_TIME)
    DURATION=$(date -u -d @${ELAPSED} +"%T")
    printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}
}
