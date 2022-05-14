#!/bin/bash

#SBATCH --job-name=05d_combine_and_genotype_gvcfs.sh
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 07:00:00
#SBATCH --mem=16000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu
#SBATCH --array=0-14

# -----------------------------------------------------------------------------
# Above array needs to contain count of elements equal to number of populations
# we're testing x number of coverage levels we're testing. e.g. If we're testing
# 3 populations and 5 coverage levels, the array should have 15 elements, 0-14.
#
# This script runs one insance for each combination of population and coverage
# level we want to analyze. For example, if we have 2 population sizes, 100 &
# 50, and three coverage levels, 50x, 30x, 15x, then this script will run
# 6 instances:
#
#   pop100_cvg50x
#   pop100_cvg30x
#   pop100_cvg15x
#   pop50_cvg50x
#   pop50_cvg30x
#   pop50_cvg15x
#
# We set the #SBATCH --array=0-5.
#
# The script needs to know which population and coverage to use in the given
# instance. We look those up in the arrays defined in init_script_vars.sh:
#
#     cvgX=(50x 30x 15x)
#     popN=(100 50)
#
# and we define cvgCnt = the number of element in the cvgX array, in this case
# 3.
#
# We map from the array index for this instance - obtained from
# $SLURM_ARRAY_TASK_ID - to the population and coverage level to use for this
# instance like so:
#
#   array    pop index =         cvg index =
#            array/cvgCnt        array%cvgCnt
#   -----    -------------       ------------
#       0          0                  0
#       1          0                  1
#       2          0                  2
#       3          1                  0
#       4          1                  1
#       5          1                  2
#
# example: for array = $SLURM_ARRAY_TASK_ID = 3,
#    the population is popN[1] = 50 and
#    the coverage is cvgX[0] = 50x
#
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05d_combine_and_genotype_gvcfs

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load gatk/4.1.9.0

# -----------------------------------------------------------------------------
# Genotype GVCFs across all samples simultaneously
# -----------------------------------------------------------------------------

# Calculate population array and coverage level array indicies for this instance
# of the array job.

cvgCnt=${#cvgX[@]}

# Set population for this instance of array job
i=$(($SLURM_ARRAY_TASK_ID / $cvgCnt))
population=${popN[i]}

# Set coverage for this instance of array job
i=$(($SLURM_ARRAY_TASK_ID % $cvgCnt))
coverage=${cvgX[i]}

# Set Map file for this sample set
MAP_FILE=${OUTPUT_DIR}/sample_pop_${population}_cvg_${coverage}_map.list

# ----------------------------------------------------------------------
# gatk Gombine GVCFS
# ----------------------------------------------------------------------

# Set Map file for this sample set
MAP_FILE=${OUTPUT_DIR}/sample_pop_${population}_cvg_${coverage}_map.list

OUT_FILE=sample_pop_${population}_cvg_${coverage}_combined_gvcfs.vcf
COMBINED_GVCFS_FILE=${OUT_FILE}

start_logging "gatk CombineGVCFs - ${OUT_FILE}"

gatk CombineGVCFs \
    -R ${REF_GENOME_FILE} \
    -V ${MAP_FILE} \
    -O ${OUTPUT_DIR}/${OUT_FILE}

stop_logging

# ----------------------------------------------------------------------
# gatk GenotypeGVCFs
# ----------------------------------------------------------------------

OUT_FILE=sample_pop_${population}_cvg_${coverage}_genotyped_gvcfs.vcf
GENOTYPED_FILE=${OUT_FILE}

start_logging "gatk GenotypeGVCFs- ${OUT_FILE}"

gatk GenotypeGVCFs \
    -R ${REF_GENOME_FILE} \
    -V ${OUTPUT_DIR}/${COMBINED_GVCFS_FILE} \
    -O ${OUTPUT_DIR}/${OUT_FILE}

stop_logging
