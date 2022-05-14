#!/bin/bash

#SBATCH --job-name=06b_run_PLINK
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu
#xSBATCH --array=0-14

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

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

# module load anaconda/3-2020.11
module load plink/1.9

# -----------------------------------------------------------------------------
# Process final filtered SNPs
# -----------------------------------------------------------------------------

# # Iterate over populations -----------------------------------------------------

# for population in ${popN[@]}; do

#     # Iterate over levels of coverage ------------------------------------------

#     for i in $(seq 0 $cvgCnt); do

# -----------------------------------------------------------------------------
# Set variables and parameters for this instance of the array job
# -----------------------------------------------------------------------------

# Set population for this instance of array job
# cvgCnt=${#cvgX[@]}
# i=$(($SLURM_ARRAY_TASK_ID / $cvgCnt))
# population=${popN[i]}
population=30

# Set coverage for this instance of array job
# i=$(($SLURM_ARRAY_TASK_ID % $cvgCnt))
# coverage=${cvgX[i]}
coverage=50x

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

# ----------------------------------------------------------------------
# Run PLINK to ID ROHs
# ----------------------------------------------------------------------

# p1=0    # Values for -homozyg-window-het
# p2=5    # Values for -homozyg-window-missing
# p3=50   # Values for -homozyg-window-snp
# p4=50   # Values for -homozyg-density
# p5=1000 # Values for -homozyg-gap
# p6=0.01 # Values for -homozyg-window-threshold
# p7=25   # Values for -homozyg-snp
# p8=10   # Values for -homozyg-kb

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

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
