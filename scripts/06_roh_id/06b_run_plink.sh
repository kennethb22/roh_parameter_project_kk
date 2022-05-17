#!/bin/bash

#SBATCH --job-name=06b_run_PLINK
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu
#SBATCH --array=0-14

# -----------------------------------------------------------------------------
# Above array needs to contain count of elements equal to number of populations
# we're testing x number of coverage levels we're testing. e.g. If we're testing
# 3 populations and 5 coverage levels, the array should have 15 elements, 0-14.
#
# This script runs one insance for each combination of population and coverage
# level we want to analyze.
#
# The script needs to know which population and coverage to use in the given
# instance. We look those values up in the arrays defined in
# init_script_vars.sh:
#
#    declare -a cvgX=(50x 30x 15x 10x 05x)
#    declare -a popN=(100 50 30)
#
# We define cvgCnt = the number of elements in the cvgX array, in this case
# 5.
#
# We map from the array index for this instance - obtained from
# $SLURM_ARRAY_TASK_ID - to the population and coverage level to use for this
# instance like so. We define the indicies to the popN and cvgX arrays as:
#
#    pop_index = array_index/cvgCnt
#    cvg_index = array_index%cvgCnt
#
# The table below shows how these functions map from the array_index to the
# actual population and coverage values contained in the cvgX and popN arrays.
#
#   array   pop     cvg     pop     cvg
#   index   index   index   value   value
#   -----   -----   -----   -----   -----
#      0       0       0     100     50x
#      1       0       1     100     30x
#      2       0       2     100     15x
#      3       0       3     100     10x
#      4       0       4     100      5x
#      5       1       0      50     50x
#      6       1       1      50     30x
#      7       1       2      50     15x
#      8       1       3      50     10x
#      9       1       4      50      5x
#     10       2       0      30     50x
#     11       2       1      30     30x
#     12       2       2      30     15x
#     13       2       3      30     10x
#     14       2       4      30      5x
#
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=06_roh_id
PREV_STEP=05_var_call
SCRIPT=06b_run_PLINK

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load plink/1.9

# -----------------------------------------------------------------------------
# Set variables and parameters for this instance of the array job
# -----------------------------------------------------------------------------

# Set population for this instance of array job
cvgCnt=${#cvgX[@]}
i=$(($SLURM_ARRAY_TASK_ID / $cvgCnt))
population=${popN[i]}

# Set coverage for this instance of array job
i=$(($SLURM_ARRAY_TASK_ID % $cvgCnt))
coverage=${cvgX[i]}

# ----------------------------------------------------------------------
# Convert VCF to PLINK format(s)
# ----------------------------------------------------------------------

IN_FILE=sample_pop_${population}_cvg_${coverage}_filtered_SNPs.vcf
OUT_FILE=sample_pop_${population}_cvg_${coverage}_plink

start_logging "Convert VCF to PLINK"

plink \
    --vcf ${INPUT_DIR}/${IN_FILE} \
    --allow-extra-chr \
    --out ${OUTPUT_DIR}/${OUT_FILE}

stop_logging

# -----------------------------------------------------------------------------
# PLINK parameters - define arrays for each of the PLINK parameters we want to
# vary and test.
# -----------------------------------------------------------------------------

declare -a phwh=(0 1 2)         # Values for -homozyg-window-het
declare -a phwm=(2 5 50)        # Values for -homozyg-window-missing
declare -a phws=(50 100 1000)   # Values for -homozyg-window-snp
declare -a phzd=(50)            # Values for -homozyg-density
declare -a phzg=(500 1000)      # Values for -homozyg-gap
declare -a phwt=(0.01 0.05 0.1) # Values for -homozyg-window-threshold
declare -a phzs=(10 100 1000)   # Values for -homozyg-snp
declare -a phzk=(100)           # Values for -homozyg-kb

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
                                IN_FILE=sample_pop_${population}_cvg_${coverage}_plink
                                OUT_FILE=sample_pop_${population}_cvg_${coverage}_plink_roh${PARAM_SUFFIX}

                                start_logging "PLINK ROH ID - ${PARAM_SUFFIX}"

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
