#!/bin/bash

# Don't need to run this one on the cluster. Just run it from the terminal.
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

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05b-2_create_subsamples

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# For each population size, create a text file with a list of individual
# sample names.
#
# 2022-03-31 - Moved creating the subsample populations to
#              01b_create_subsample_lists
# -----------------------------------------------------------------------------

# for population in ${popN[@]}; do

#     # Randomly choose n individuals from original population. Save list of
#     # individuals to a text file
#     echo ${population}
#     echo ${SAMPLE_ID_LIST}
#     echo ${OUTPUT_DIR}/sample_id_list_pop_${population}.txt

#     shuf -n ${population} ${SAMPLE_ID_LIST} >${OUTPUT_DIR}/sample_id_list_pop_${population}.txt

# done

# -----------------------------------------------------------------------------
# For each population size and coverage level, create a map file with a list
# of full paths to each file in that population and coverage level group.
# -----------------------------------------------------------------------------
# Iterate over populations -----------------------------------------------------

declare -a popN=(100 50 30)
declare -a cvgX=(50x 30x 15x 10x 05x)

cvgCnt=${#cvgX[@]}
let cvgCnt-=1
for population in ${popN[@]}; do

    # Iterate over levels of coverage ------------------------------------------

    for i in $(seq 0 $cvgCnt); do

        # Set map file

        MAP_FILE=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}_map.list
        if [ -f "$MAP_FILE" ]; then
            rm ${MAP_FILE}
        fi

        # Set path to sample vcf files created by gatk haplotypeCaller

        CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_cvg_${cvgX[i]}

        # Create map file for this level op
        while read -a line; do

            VCF_FILE=${line[0]}_cvg_${cvgX[i]}.g.vcf

            echo ${CVG_OUTPUT_DIR}/${VCF_FILE} >>${MAP_FILE}

        done <${OUTPUT_DIR}/sample_id_list_pop_${population}.txt
    done
done
