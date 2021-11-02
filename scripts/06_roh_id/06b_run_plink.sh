#!/bin/bash
#
#  +-----------------+
#  | REQUIRES 20 CPU |
#  +-----------------+
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
# Set Variables
# -----------------------------------------------------------------------------

## Set step name

STEP=06_roh_id
PREV_STEP=05_var_call
SCRIPT=06b_run_PLINK

## Load variables and functions from settings file

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
## Load modules
# -----------------------------------------------------------------------------

module load anaconda/3-2020.11

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
# Process final filtered SNPs
# -----------------------------------------------------------------------------

# Iterate over populations -----------------------------------------------------

for population in ${popN[@]}; do

    # Iterate over levels of coverage ------------------------------------------

    for i in $(seq 0 $cvgCnt); do

        # ----------------------------------------------------------------------
        # Convert VCF to PLINK format(s)
        # ----------------------------------------------------------------------

        ## Set action and start time for log

        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_final_SNPs.vcf
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_plink

        start_logging "Convert VCF to PLINK"

        plink \
            --vcf ${INPUT_DIR}/${IN_FILE} \
            --allow-extra-chr \
            --out ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging

        # ----------------------------------------------------------------------
        # Run PLINK to ID ROHs
        # ----------------------------------------------------------------------

        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_plink
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_plink_roh

        start_logging "PLINK ROH ID"

        plink --bim ${OUTPUT_DIR}/${IN_FILE}.bim \
            --bed ${OUTPUT_DIR}/${IN_FILE}.bed \
            --fam ${OUTPUT_DIR}/${IN_FILE}.fam \
            --allow-extra-chr --no-pheno \
            --homozyg-window-het 1 \
            --homozyg-window-missing 5 \
            --homozyg-window-snp 50 \
            --homozyg-density 50 \
            --homozyg-gap 1000 \
            --homozyg-window-threshold 0.05 \
            --homozyg-snp 25 \
            --homozyg-kb 10 \
            --out ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging
    done

done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
