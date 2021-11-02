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
# Define Functions
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set Variables
# -----------------------------------------------------------------------------

## Load variables from settings file
source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

## Set step name
STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05b_genomicsdb

## Create variable pointing to user's home directory for this project and step
HOME_STEP_DIR=/home/${USER}/${PROJECT}/data/${STEP}/

## Set the input directory, which is the output directory from the previous
## step
INPUT_DIR=/scratch/${USER}/${PROJECT}/data/${STEP}/output

## Create a working directory on /scratch
mkdir -p /scratch/${USER}/${PROJECT}/data/${STEP}
chmod -R 700 /scratch/${USER}
WORK_DIR=/scratch/${USER}/${PROJECT}/data/${STEP}

## Set output directory
OUTPUT_DIR=${WORK_DIR}/output

# -----------------------------------------------------------------------------
# Load Modules
# -----------------------------------------------------------------------------
module load gatk/4.1.4.0
# module load picard/1.79
# module load vcftools/0.1.14

# -----------------------------------------------------------------------------
# Create log file
# -----------------------------------------------------------------------------

# Set name of log file to track execution times
LOG_FILE=${OUTPUT_DIR}/${SCRIPT}_log.txt

# Delete log file if it exists

if [ -f "$LOG_FILE" ]; then
    rm ${LOG_FILE}
fi

# Write header to log file

printf "%-80s   %8s   %8s   %8s\n" "Action - Output" "Start" "End" "Duration" >>${LOG_FILE}

# -----------------------------------------------------------------------------
# Genotype GVCFs across all samples simultaneously
# -----------------------------------------------------------------------------

# Iterate over each combination of population and coverage level. Sets of
# population and coverage are set in init_script_vars.sh

# Iterate over populations -----------------------------------------------------

for population in ${popN[@]}; do

    # Iterate over levels of coverage ------------------------------------------

    for i in $(seq 0 $cvgCnt); do

        ## Set Map file for this sample set
        MAP_FILE=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}_map.list

        # ----------------------------------------------------------------------
        # gatk Gombine GVCFS
        # ----------------------------------------------------------------------

        ## Set action and start time for log

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_combined_gcvfs.vcf
        COMBINED_GVCFS_FILE=${OUT_FILE}
        ACTION="gatk CombineGVCFS - ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        ## Do GenomicsDBImport

        gatk CombineGVCFs \
            -R ${REF_GENOME_FILE} \
            -V ${MAP_FILE} \
            -O ${OUTPUT_DIR}/${OUT_FILE}

        ## Write line to log file for this action

        END_TIME_HHMM=$(date +"%T")
        END_TIME=$(date +%s)
        ELAPSED=$(expr $END_TIME - $START_TIME)
        DURATION=$(date -u -d @${ELAPSED} +"%T")
        printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}

        # ----------------------------------------------------------------------
        # gatk GenotypeGVCFs
        # ----------------------------------------------------------------------

        ## Set action and start time for log

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_genotyped_gcvfs.vcf
        GENOTYPED_FILE=${OUT_FILE}
        ACTION="gatk GenotypeGCVFS - ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        ## Do GenotypeGCVFs

        gatk GenotypeGVCFs \
            -R ${REF_GENOME_FILE} \
            -V ${OUTPUT_DIR}/${COMBINED_GVCFS_FILE} \
            -O ${OUTPUT_DIR}/${OUT_FILE}

        ## Write line to log file for this action

        END_TIME_HHMM=$(date +"%T")
        END_TIME=$(date +%s)
        ELAPSED=$(expr $END_TIME - $START_TIME)
        DURATION=$(date -u -d @${ELAPSED} +"%T")
        printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}

        # ----------------------------------------------------------------------
        # Flag Variants that do not pass filters
        # ----------------------------------------------------------------------

        ## Set action and start time for log

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_gcvfs.vcf
        FILTERED_FILE=${OUT_FILE}
        ACTION="gatk VariantFiltration- ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        ## Do VariantFiltration

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

        ## Write line to log file for this action

        END_TIME_HHMM=$(date +"%T")
        END_TIME=$(date +%s)
        ELAPSED=$(expr $END_TIME - $START_TIME)
        DURATION=$(date -u -d @${ELAPSED} +"%T")
        printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}

        # ----------------------------------------------------------------------
        # Keep only SNPs that pass the above filters
        # ----------------------------------------------------------------------

        ## Set action and start time for log

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_final_SNPs.vcf
        ACTION="gatk SelectVariants - ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        # Do SelectVariants

        # Keep only SNPs that pass the above filters
        gatk SelectVariants \
            -R ${REF_GENOME_FILE} \
            -V ${OUTPUT_DIR}/${FILTERED_FILE} \
            --select-type-to-include SNP \
            -select 'vc.isNotFiltered()' \
            -O ${OUTPUT_DIR}/${OUT_FILE}

        ## Write line to log file for this action

        END_TIME_HHMM=$(date +"%T")
        END_TIME=$(date +%s)
        ELAPSED=$(expr $END_TIME - $START_TIME)
        DURATION=$(date -u -d @${ELAPSED} +"%T")
        printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}

        ## Delete GenomicsDB for this sample group

    done

done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
