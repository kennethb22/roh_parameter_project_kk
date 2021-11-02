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
STEP=01_slim

# -----------------------------------------------------------------------------
# These arrays are used in the steps:
#    04_downsample.sh
#    05a_snp_calling.sh
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
# declare -a popN=(03)
declare -a popN=(01)

# -----------------------------------------------------------------------------
## Set location of reference genome file
REF_GENOME_FILE_PATH=/scratch/${USER}/${PROJECT}/data/${STEP}/output/slim_output_files_m5e-07_r1e-8_p500

REF_GENOME_FILE_NAME=ancestral

REF_GENOME_FILE=${REF_GENOME_FILE_PATH}/${REF_GENOME_FILE_NAME}.fasta
