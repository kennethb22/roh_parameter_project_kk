#!/bin/bash
#
#  +----------------------+
#  | REQUEST 4 CPU + 60gb |
#  +----------------------+
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
SCRIPT=06a_run_bcftoolsROH

## Load variables and functions from settings file

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
## Load modules
# -----------------------------------------------------------------------------

module load bcftools/1.10.2
module load samtools/1.11

## Create arrays of downsample levels to be used.
## Coverage level in NNx for display in output file names
declare -a cvgX=(50x 30x 15x 10x 05x)

## Coverage level fraction to supply to samtools
declare -a cvgP=(1.0 0.6 0.3 0.2 0.1)

## Get length of the coverage level arrays. Subtract 1 because arrays are zero
## based, and we'll iterate over the arrays from 0 to cvgCnt
cvgCnt=${#cvgX[@]}
let cvgCnt-=1

## Create array of population sizes we want to test
declare -a popN=(100 50 30)

# -----------------------------------------------------------------------------
# Process final filtered SNPs
# -----------------------------------------------------------------------------

# Iterate over populations -----------------------------------------------------

for population in ${popN[@]}; do

    # Iterate over levels of coverage ------------------------------------------

    for i in $(seq 0 $cvgCnt); do

        # ----------------------------------------------------------------------
        # Generate allele frequency files for use with bcftools roh & index
        # with tabix
        # ----------------------------------------------------------------------

        ## Set action and start time for log

        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_SNPs.vcf
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_filteredSNPs.tab.gz

        start_logging "bcftools query/tabix  - ${OUT_FILE}"

        # Extract allele frequencies from sample vcf

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
            ${INPUT_DIR}/${IN_FILE} | bgzip -c >${OUTPUT_DIR}/${OUT_FILE}

        tabix -s1 -b2 -e2 ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging

        # ----------------------------------------------------------------------
        # ROH Calling - Using Genotypes only
        # ----------------------------------------------------------------------

        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_SNPs.vcf
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_roh_gt.txt
        AF_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_filteredSNPs.tab.gz

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

        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_roh_gt.txt
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_roh_gt_RG_ONLY.txt
        grep "^RG" ${OUTPUT_DIR}/${IN_FILE} > \
            ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging

        # ----------------------------------------------------------------------
        # ROH Calling - Using Genotypes and genotype liklihoods
        # ----------------------------------------------------------------------

        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_SNPs.vcf
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_roh_pl.txt
        AF_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_filteredSNPs.tab.gz

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

        IN_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_roh_pl.txt
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_bcftools_roh_pl_RG_ONLY.txt
        grep "^RG" ${OUTPUT_DIR}/${IN_FILE} > \
            ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging

    done

done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
