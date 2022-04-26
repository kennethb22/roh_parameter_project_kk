#!/bin/bash
#
#  This script takes the list of individual IDs that were created by SLiM and
#  uses that list to create randomly selected subpopulations, using the
#  population sizes specified in init_script_vars.sh. It first selects the
#  individuals for the largest population, then selects the individuals for the
#  second largest population from that list, then pulls the individuals for the
#  third largest population from the second largest, etc.
#
#  This script does not need to be submitted to the queue. Just run manaully in
#  the shell.
#
#  This script must be made executable like this
#    chmod +x my_script
#
#  VERY IMPORTANT: Only run this script ONCE after you've created the sample
#  populations in SLiM. It will overwrite any existing subpopulation lists.
#
#  TO DO: add a check to see if subpopulation lists with the given names
#         already exist, stop if they do, and alert the user.

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=01b_create_subsample_lists
PREV_STEP=01_run_slim
SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------

for population in ${popN[@]}; do

    # Randomly choose n individuals from original population. Save list of
    # individuals to a text file
    # echo PREV: ${PREV_POPULATION}
    # echo POP : ${population}
    # echo POPN: ${popN[0]}
    # echo cmp: ${population} == ${popN[0]}
    # echo INFILE: ${SAMPLE_ID_LIST}
    # echo ${INIT_OUTPUT_DIR}/sample_id_list_pop_${population}.txt

    if [ "$population" -eq "${popN[0]}" ]; then
        # echo "First element"
        SOURCE_LIST=${BASE_SAMPLE_ID_LIST}
    else
        # echo "Not first element"
        SOURCE_LIST=${INIT_OUTPUT_DIR}/sample_id_list_pop_${PREV_POPULATION}.txt
    fi

    # echo "Source List: "${SOURCE_LIST}

    PREV_POPULATION=${population}
    shuf -n ${population} ${SOURCE_LIST} >${INIT_OUTPUT_DIR}/sample_id_list_pop_${population}.txt

done
