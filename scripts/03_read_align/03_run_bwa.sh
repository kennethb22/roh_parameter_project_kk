#!/bin/bash

#SBATCH --job-name=03_read_align
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH --mem=20000
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=kbk0024@auburn.edu

# srun -N1 -n20 -t02:00:00 --mem=20000 --pty bash

#   +-----------------------+
#   |  USE:                 |
#   |    - LARGE queue      |
#   |    - 20 CPU + 20 Gb   |
#   +-----------------------+
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=03_read_align
PREV_STEP=02_read_sim_and_qc
SCRIPT=03_run_bwa.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /scratch/kbk0024/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load bwa/0.7.17
module load samtools/1.11

# -----------------------------------------------------------------------------
# Index reference genomes
# -----------------------------------------------------------------------------

# TO DO: move this indexing to step O1  - make it step 01c_index_ref_genome

# bwa index ${REF_GENOME_FILE}

# -----------------------------------------------------------------------------
# Align reads to reference genome
# -----------------------------------------------------------------------------

cd ${OUTPUT_DIR}

while read -a line; do

    OUT_FILE=${line[0]}_genome.bam
    start_logging "bwa align - ${OUT_FILE}"

    bwa mem -t 20 -M \
        ${REF_GENOME_FILE} \
        ${INPUT_DIR}/${line[0]}_f.fq ${INPUT_DIR}/${line[0]}_r.fq >${OUT_FILE}

    stop_logging

done <${SAMPLE_ID_LIST}

# -----------------------------------------------------------------------------
# Sort aligned read bam files
# -----------------------------------------------------------------------------

while read -a line; do

    OUT_FILE=${line[0]}_genome_sorted.bam
    start_logging "samtools sort - ${OUT_FILE}"

    samtools sort -@ 19 \
        -o ${OUT_FILE} \
        ${line[0]}_genome.bam

    stop_logging

done <${SAMPLE_ID_LIST}

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
