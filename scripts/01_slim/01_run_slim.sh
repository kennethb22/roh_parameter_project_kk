#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - MEDIUM queue     |
#   |    - 2 CPU + 6 Gb     |
#   +-----------------------+
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
#  My preferred setup before running:
#    -- script to be run in /home/scripts
#    -- project directory (of same name as script) in /home/
#    -- /input/ and /output/ subdirs within project dir
#
# IMPORTANT: before running this script, you must set SLiM parameters in
# init_script_vars.sh
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=01_slim
SCRIPT=01_run_slim

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load anaconda/3-2020.11

# -----------------------------------------------------------------------------
# Run SLiM - SLIM_OUT_DIR, parameters, and SLIM_PARAM_FILE are set in
# init_script_vars.sh
# -----------------------------------------------------------------------------

mkdir ${OUTPUT_DIR}/${SLIM_OUT_DIR}

start_logging "Run SLiM - ${SLIM_OUT_DIR}"

# models chromosome with coding and non-coding regions
slim \
    -d POP_SIZE=${POP_SIZE} \
    -d MUTATION_RATE=${MUTATION_RATE} \
    -d RECOMB_RATE=${RECOMB_RATE} \
    -d "OUT_PATH='${OUTPUT_DIR}/${SLIM_OUT_DIR}/'" \
    ${SLIM_PARAM_FILE}

stop_logging

# mail -s 'SLiM run finished - Starting BCFtools' ${EMAIL} <<<'SLiM run finished - Starting BCFtools'

# -----------------------------------------------------------------------------
# Format SLiM output for read simulation - FILE_LABELS and SAMPLE_ID_LIST and
# FASTA_OUT_DIR are set in init_script_vars.sh. FASTA_OUT_DIR is used as the
# input directory in 02a_run_art.sh
# -----------------------------------------------------------------------------

module purge
module load bcftools/1.13

start_logging "Format SLiM Output for read simulation - ${SLIM_OUT_DIR}"

VCF_OUT_DIR=${OUTPUT_DIR}/sample_vcf_files_${FILE_LABELS}
mkdir ${VCF_OUT_DIR}
VCF_FILE_LIST=${OUTPUT_DIR}/vcf_file_list_${FILE_LABELS}.txt

## Compress output VCF
bgzip ${OUTPUT_DIR}/${SLIM_OUT_DIR}/final_pop.vcf

## Index compressed VCF
tabix -f ${OUTPUT_DIR}/${SLIM_OUT_DIR}/final_pop.vcf.gz

## Split output VCF into sample-specific VCF files

bcftools +split -O z -o ${VCF_OUT_DIR}/ ${OUTPUT_DIR}/${SLIM_OUT_DIR}/final_pop.vcf.gz

## Randomly select sample VCF files for conversation to FASTAs
## >> Selecting 100, can be downsampled later

ls ${VCF_OUT_DIR} i*.vcf.gz | sort -R | tail -100 >${VCF_FILE_LIST}
sed 's/.vcf.gz//g' ${VCF_FILE_LIST} >${SAMPLE_ID_LIST}

## Convert sample VCF files to two separate haplotype FASTAs per individual

while read -a line; do
    tabix ${VCF_OUT_DIR}/${line[0]}.vcf.gz

    bcftools norm --check-ref s \
        --fasta-ref ${REF_GENOME_FILE} \
        --multiallelics - \
        --do-not-normalize \
        --output ${VCF_OUT_DIR}/norm_${line[0]}.vcf \
        /${VCF_OUT_DIR}/${line[0]}.vcf.gz

    bgzip ${VCF_OUT_DIR}/norm_${line[0]}.vcf
    tabix ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

    bcftools +fixref ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz -- \
        --fasta-ref ${REF_GENOME_FILE}

    bcftools consensus --haplotype 1 \
        --fasta-ref ${REF_GENOME_FILE} \
        --output ${FASTA_OUT_DIR}/${line[0]}_1.fasta \
        ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

    bcftools consensus --haplotype 2 \
        --fasta-ref ${REF_GENOME_FILE} \
        --output ${FASTA_OUT_DIR}/${line[0]}_2.fasta \
        ${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

done <${SAMPLE_ID_LIST}

stop_logging

# mail -s 'SLiM run finished - submit Rviz' ${EMAIL} <<<'SLiM run finished'

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
