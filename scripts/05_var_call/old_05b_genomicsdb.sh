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
#    run_script_scratch my_script.sh
#

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05b_genomicsdb

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------

module load gatk/4.1.4.0

# -----------------------------------------------------------------------------
# Genotype GVCFs across all samples simultaneously
# -----------------------------------------------------------------------------

# Iterate over each combination of population and coverage level. Sets of
# population and coverage are set in init_script_vars.sh

# Iterate over populations -----------------------------------------------------

for population in ${popN[@]}; do

    # Iterate over levels of coverage ------------------------------------------

    for i in $(seq 0 $cvgCnt); do

        # Set Map file for this sample set
        MAP_FILE=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}_map.list

        # ----------------------------------------------------------------------
        # gatk Gombine GVCFS
        # ----------------------------------------------------------------------

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_combined_gcvfs.vcf
        COMBINED_GVCFS_FILE=${OUT_FILE}

        start_logging "gatk CombineGVCFS - ${OUT_FILE}"

        # Do GenomicsDBImport

        gatk CombineGVCFs \
            -R ${REF_GENOME_FILE} \
            -V ${MAP_FILE} \
            -O ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging

        # ----------------------------------------------------------------------
        # gatk GenotypeGVCFs
        # ----------------------------------------------------------------------

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_genotyped_gcvfs.vcf
        GENOTYPED_FILE=${OUT_FILE}

        start_logging "gatk GenotypeGCVFS - ${OUT_FILE}"

        # Do GenotypeGCVFs

        gatk GenotypeGVCFs \
            -R ${REF_GENOME_FILE} \
            -V ${OUTPUT_DIR}/${COMBINED_GVCFS_FILE} \
            -O ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging

        # ----------------------------------------------------------------------
        # Flag Variants that do not pass filters
        # ----------------------------------------------------------------------

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_gcvfs.vcf
        FILTERED_FILE=${OUT_FILE}

        start_logging "gatk VariantFiltration- ${OUT_FILE}"

        # Do VariantFiltration

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

        # ----------------------------------------------------------------------
        # Keep only SNPs that pass the above filters
        # ----------------------------------------------------------------------

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_final_SNPs.vcf

        start_logging "gatk SelectVariants - ${OUT_FILE}"

        # Do Select Variants

        gatk SelectVariants \
            -R ${REF_GENOME_FILE} \
            -V ${OUTPUT_DIR}/${FILTERED_FILE} \
            --select-type-to-include SNP \
            -select 'vc.isNotFiltered()' \
            -O ${OUTPUT_DIR}/${OUT_FILE}

        stop_logging

    done

done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
