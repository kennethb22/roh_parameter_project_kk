#!/bin/bash

#SBATCH --job-name=05e_filter_variants
#SBATCH -N 1
#SBATCH -n 2
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

#
#  NOTES - Run on Alabama Supercomputer Cluster:
#
#  large queue:  REQUEST 2 CPU + 60 gb
#
#  2021-11-14 - 05d_genotype_and_filter_p100_10x completed but did't
#    produce filtered.vcf output file. Checked job output log, execution just
#    stopped during the CombineVCFS step. I think that happened because I had
#    exceeded my disk quota. Scripts were making a new backup of the output
#    directory from scratch for every copy of 05d being run. Filled up quota
#    quickly. Disabling backup code in the script. Will manually backup of
#    output directory from scratch to home after all 05d scripts finish.
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05e_filter_variants

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load gatk/4.1.9.0

# -----------------------------------------------------------------------------
# Set variables and parameters for this instance of the array job
# -----------------------------------------------------------------------------

# ## Create arrays of downsample levels to be used.
# ## Coverage level in NNx for display in output file names
# declare -a cvgX=(50x 30x 15x 10x 05x)

# ## Coverage level fraction to supply to samtools
# declare -a cvgP=(1.0 0.6 0.3 0.2 0.1)

# ## Get length of the coverage level arrays. Subtract 1 because arrays are zero
# ## based, and we'll iterate over the arrays from 0 to cvgCnt
# cvgCnt=${#cvgX[@]}
# let cvgCnt-=1

# ## Create array of population sizes we want to test
# # declare -a popN=(100 50 30)
# declare -a popN=(100 50 30)

# Calculate population array and coverage level array indicies for this instance
# of the array job.

cvgCnt=${#cvgX[@]}

# Set population for this instance of array job
i=$(($SLURM_ARRAY_TASK_ID / $cvgCnt))
population=${popN[i]}

# Set coverage for this instance of array job
i=$(($SLURM_ARRAY_TASK_ID % $cvgCnt))
coverage=${cvgX[i]}

# ------------------------------------------------------------------------------
# Flag Variants that do not pass filters
# ------------------------------------------------------------------------------

OUT_FILE=sample_pop_${population}_cvg_${coverage}_filtered_gvcfs.vcf
FILTERED_FILE=${OUT_FILE}
GENOTYPED_FILE=sample_pop_${population}_cvg_${coverage}_genotyped_gvcfs.vcf

start_logging "gatk VariantFiltration- ${OUT_FILE}"

gatk VariantFiltration \
    -R ${REF_GENOME_FILE} \
    -V ${OUTPUT_DIR}/${GENOTYPED_FILE} \
    -O ${OUTPUT_DIR}/${OUT_FILE} \
    --filter-name "QD" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS" \
    --filter-expression "FS > 40.0" \
    --filter-name "SOR" \
    --filter-expression "SOR > 5.0" \
    --filter-name "MQ" \
    --filter-expression "MQ < 20.0" \
    --filter-name "MQRankSum" \
    --filter-expression " MQRankSum < -3.0 || MQRankSum > 3.0" \
    --filter-name "ReadPosRankSum" \
    --filter-expression "ReadPosRankSum < -8.0"

stop_logging

# ------------------------------------------------------------------------------
# Keep only SNPs that pass the above filters
# ------------------------------------------------------------------------------

OUT_FILE=sample_pop_${population}_cvg_${coverage}_filtered_SNPs.vcf

start_logging "gatk SelectVariants - ${OUT_FILE}"

gatk SelectVariants \
    -R ${REF_GENOME_FILE} \
    -V ${OUTPUT_DIR}/${FILTERED_FILE} \
    --select-type-to-include SNP \
    -select 'vc.isNotFiltered()' \
    -O ${OUTPUT_DIR}/${OUT_FILE}

stop_logging

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
