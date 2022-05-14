#!/bin/bash

#SBATCH --job-name=06a_run_bcftoolsROH
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH --mem=32000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu
#SBATCH --array=0-14

#  NOTES: Running on Alabama Supercomputer Cluster
#
#  +----------------------+
#  | REQUEST 4 CPU + 60gb |
#  +----------------------+
#

# -----------------------------------------------------------------------------
# Set Variables
# -----------------------------------------------------------------------------

STEP=06_roh_id
PREV_STEP=05_var_call
SCRIPT=06a_run_bcftoolsROH

## Load variables and functions from settings file

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
## Load modules
# -----------------------------------------------------------------------------

module load bcftools/1.11
module load samtools/1.11

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

# -----------------------------------------------------------------------------
# Process final filtered SNPs
# -----------------------------------------------------------------------------

# ----------------------------------------------------------------------
# Generate allele frequency files for use with bcftools roh & index
# with tabix
# ----------------------------------------------------------------------

IN_FILE=sample_pop_${population}_cvg_${coverage}_filtered_SNPs.vcf
OUT_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_filteredSNPs.tab.gz

start_logging "bcftools query/tabix  - ${OUT_FILE}"

# Extract allele frequencies from sample vcf

bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
    ${INPUT_DIR}/${IN_FILE} | bgzip -c >${OUTPUT_DIR}/${OUT_FILE}

tabix -s1 -b2 -e2 ${OUTPUT_DIR}/${OUT_FILE}

stop_logging

# ----------------------------------------------------------------------
# ROH Calling - Using Genotypes only
# ----------------------------------------------------------------------

IN_FILE=sample_pop_${population}_cvg_${coverage}_filtered_SNPs.vcf
OUT_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_roh_gt.txt
AF_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_filteredSNPs.tab.gz

start_logging "bcftools roh gt - ${OUT_FILE}"

bcftools roh \
    --GTs-only 30 \
    --threads 20 \
    -o ${OUTPUT_DIR}/${OUT_FILE} \
    --AF-file ${OUTPUT_DIR}/${AF_FILE} \
    ${INPUT_DIR}/${IN_FILE}

# ----------------------------------------------------------------------
# Extract information on ROHS i.e., exlcude information on individual
# sites
# ----------------------------------------------------------------------

IN_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_roh_gt.txt
OUT_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_roh_gt_RG_ONLY.txt
grep "^RG" ${OUTPUT_DIR}/${IN_FILE} > \
    ${OUTPUT_DIR}/${OUT_FILE}

stop_logging

# ----------------------------------------------------------------------
# ROH Calling - Using Genotypes and genotype liklihoods
# ----------------------------------------------------------------------

IN_FILE=sample_pop_${population}_cvg_${coverage}_filtered_SNPs.vcf
OUT_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_roh_pl.txt
AF_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_filteredSNPs.tab.gz

start_logging "bcftools roh pl - ${OUT_FILE}"

bcftools roh \
    --threads 20 \
    -o ${OUTPUT_DIR}/${OUT_FILE} \
    --AF-file ${OUTPUT_DIR}/${AF_FILE} \
    ${INPUT_DIR}/${IN_FILE}

# ----------------------------------------------------------------------
# Extract information on ROHS i.e., exlcude information on individual
# sites
# ----------------------------------------------------------------------

IN_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_roh_pl.txt
OUT_FILE=sample_pop_${population}_cvg_${coverage}_bcftools_roh_pl_RG_ONLY.txt
grep "^RG" ${OUTPUT_DIR}/${IN_FILE} > \
    ${OUTPUT_DIR}/${OUT_FILE}

stop_logging

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
