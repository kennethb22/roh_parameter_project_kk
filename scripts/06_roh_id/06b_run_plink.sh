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
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=06_roh_id
PREV_STEP=05_var_call
SCRIPT=06b_run_PLINK

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load anaconda/3-2020.11

## Create arrays of downsample levels to be used.
## Coverage level in NNx for display in output file names
# declare -a cvgX=(50x 30x 15x 10x 05x)
# declare -a cvgX=(50x 30x)
declare -a cvgX=(30x)

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
declare -a popN=(100 50 30)
# declare -a popN=(01)

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
        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_SNPs.vcf
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_plink

        start_logging "Convert VCF to PLINK"

        plink \
            --vcf ${INPUT_DIR}/${IN_FILE} \
            --allow-extra-chr \
            --out ${OUTPUT_DIR}/${OUT_FILE}

        # ----------------------------------------------------------------------
        # Run PLINK to ID ROHs
        # ----------------------------------------------------------------------

        for p1 in ${phwh[@]}; do                             # -homozyg-window-het
            for p2 in ${phwm[@]}; do                         # -homozyg-winodow-missing
                for p3 in ${phws[@]}; do                     # -homozyg-window-snp
                    for p4 in ${phzd[@]}; do                 # -homozyg-density
                        for p5 in ${phzg[@]}; do             # -homozyg-gap
                            for p6 in ${phwt[@]}; do         # -homozyg-window threshold
                                for p7 in ${phzs[@]}; do     # -homozyg-snp
                                    for p8 in ${phzk[@]}; do # -homozyg-kb

                                        PARAM_SUFFIX=_phwh_${p1}_phwm_${p2}_phws_${p3}_phzd_${p4}_phzg_${p5}_phwt_${p6}_phzs_${p7}_phzk_${p8}
                                        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_plink
                                        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_plink_roh${PARAM_SUFFIX}

                                        start_logging "PLINK ROH ID"

                                        plink --bim ${OUTPUT_DIR}/${IN_FILE}.bim \
                                            --bed ${OUTPUT_DIR}/${IN_FILE}.bed \
                                            --fam ${OUTPUT_DIR}/${IN_FILE}.fam \
                                            --homozyg-window-het ${p1} \
                                            --homozyg-window-missing ${p2} \
                                            --homozyg-window-snp ${p3} \
                                            --homozyg-density ${p4} \
                                            --homozyg-gap ${p5} \
                                            --homozyg-window-threshold ${p6} \
                                            --homozyg-snp ${p7} \
                                            --homozyg-kb ${p8} \
                                            --out ${OUTPUT_DIR}/${OUT_FILE}

                                        stop_logging
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done

    done
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
