#!/bin/bash
#
#SBATCH --job-name=05a_run_downsample
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 03:30:00
#SBATCH --mem=20000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu
#SBATCH --array=0-4

# IMPORTANT: The SBATCH --array parameter above needs to match the number of
# elements in the levels of coverage array cvgX set in init_scripts_vars.sh. For
# instance, if cvgX has 5 coverage levels, set --array=0-4. If cvgX has 4
# coverage levels, set --array=0-3.
#
# --array=0-4 tells Slrum to create 5 instances of this job. Each job will have
# a different array index from the set 0, 1, 2, 3, 4. The array index for each
# instance of the job is stored in $SLURM_ARRAY_TASK_ID. In this script,
# $SLURM_ARRAY_TASK_ID is used to look up the coverage level in cvgX in order
# to create the input and output directory names for each job. This script
# creates a separate sbatch job for each coverage level in cvgX.
#
# 2022-04-28 - Modify script to run as an sbatch job array on Easley.
#            - Update calls to picard to use new command line syntax.
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
# SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")
SCRIPT=05a_add_read_groups.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
# sleep 5

module load picard/2.23.9
# module load gatk/4.1.4.0

# -----------------------------------------------------------------------------
# Run Haplotype caller on all sample files
# -----------------------------------------------------------------------------

# Create output directory for each coverage level.

CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}
mkdir ${CVG_OUTPUT_DIR}
echo CVG_OUTPUT_DIR: ${CVG_OUTPUT_DIR}

# Set input directory for this coverage

CVG_INPUT_DIR=${INPUT_DIR}/sample_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}
echo CVG_INPUT_DIR: ${CVG_INPUT_DIR}

# Process each individual sample file ----------------------------------

while read -a line; do

    # # ------------------------------------------------------------------
    # # Add read group information
    # # ------------------------------------------------------------------

    IN_FILE=${line[0]}_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}.bam
    OUT_FILE=${line[0]}_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}_rgroups.bam
    RGROUPS_FILE=${OUT_FILE}
    echo IN_FILE: ${IN_FILE}
    echo OUT_FILE: ${OUT_FILE}
    start_logging "picard/AddOrReplaceReadGroups - ${OUT_FILE}"

    ## Do picard read groups

    # Format for Alabama Supercomputer
    # java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar \

    # Format for Easley.auburn.edu
    java -jar /tools/picard-2.23.9/libs/picard.jar AddOrReplaceReadGroups \
        I=${CVG_INPUT_DIR}/${IN_FILE} \
        O=${CVG_OUTPUT_DIR}/${OUT_FILE} SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=${line[0]} RGPU=unit4

    stop_logging

    # ------------------------------------------------------------------
    # Index BAM files
    # ------------------------------------------------------------------

    OUT_FILE=${line[0]}_cvg_${cvgX[$SLURM_ARRAY_TASK_ID]}_rgroups.bai
    echo BAM_OUT_FILE: ${OUT_FILE}
    start_logging "picard/BuildBamIndex - ${OUT_FILE}"

    ## Run picard BuildBamIndex

    # Format for Alabama Supercomputer
    # java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/BuildBamIndex.jar \

    # Format for Easley.auburn.edu
    java -jar /tools/picard-2.23.9/libs/picard.jar BuildBamIndex \
        I=${CVG_OUTPUT_DIR}/${RGROUPS_FILE}

    stop_logging

done <${SAMPLE_ID_LIST}

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
