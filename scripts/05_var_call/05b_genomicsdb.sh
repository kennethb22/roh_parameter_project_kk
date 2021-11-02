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
        MAP_FILE=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}_map.txt

        ## Genomics DBImport

        ## Set action and start time for log

        OUT_FILE=pop_${population}_cvg_${cvgX[i]}
        DB_FILE=${OUTPUT_DIR}/gatk_genomicsdb_${OUT_FILE}
        ACTION="gatk GenomicsDBImport - ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        ## Do GenomicsDBImport

        # gatk GenomicsDBImport \
        #     -L 1 \
        #     --sample-name-map ${MAP_FILE} \
        #     --genomicsdb-workspace-path ${DB_FILE}

        OUT_DIR=/scratch/aubkbk001/roh_param_project/data/05_var_call/output/sample_pop_03_cvg_05x

        gatk GenomicsDBImport \
            -V ${OUT_DIR}/i31_pop_03_cvg_05x.g.vcf \
            -V ${OUT_DIR}/i322_pop_03_cvg_05x.g.vcf \
            -V ${OUT_DIR}/i428_pop_03_cvg_05x.g.vcf \
            -L 1 \
            --genomicsdb-workspace-path simchr_database_gatk_genomics

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

        OUT_FILE=pop_${population}_cvg_${cvgX[i]}
        ACTION="gatk GenotypeGCVFS - ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        ## Do GenotypeGCVFs

        sleep 1

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

        OUT_FILE=pop_${population}_cvg_${cvgX[i]}
        ACTION="gatk VariantFiltration- ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        ## Do GenotypeGCVFs

        sleep 1

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

        OUT_FILE=${line[0]}_pop_${population}_cvg_${cvgX[i]}
        ACTION="gatk SelectVariants - ${OUT_FILE}"
        START_TIME_HHMM=$(date +" %T")
        START_TIME=$(date +%s)

        ## Do GenotypeGCVFs

        sleep 1

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

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh

# ## --------------------------------
# ## Genotype GVCFs across all samples simultaneously

# ## -V [sample VCF file]; unfortunately, no easy way to just provide a file with names;
# ## will need one -V line per sample
# gatk GenomicsDBImport \
#     -V sampfix_3811_S67_nobaseQrecal.g.vcf \
#     -V sampfix_3850_S45_nobaseQrecal.g.vcf \
#     -V sampfix_3910_S63_nobaseQrecal.g.vcf \
#     --genomicsdb-workspace-path simchr_database_gatk_genomics

# gatk GenotypeGVCFs \
#     -R $ref \
#     -V gendb://simchr_database_gatk_genomics \
#     -O simchr_genotype_output_nobaseQrecal.vcf

# ## Flag variants that do not pass filters
# gatk VariantFiltration \
#     -R $ref \
#     -V simchr_genotype_output_nobaseQrecal.vcf \
#     -O simchr_genotype_output_nobaseQrecal_filtered.vcf \
#     --filter-name "QD" \
#     --filter-expression "QD < 2.0" \
#     --filter-name "FS" \
#     --filter-expression "FS > 40.0" \
#     --filter-name "SOR" \
#     --filter-expression "SOR > 5.0" \
#     --filter-name "MQ" \
#     --filter-expression "MQ < 20.0" \
#     --filter-name "MQRankSum" \
#     --filter-expression " MQRankSum < -3.0 || MQRankSum > 3.0" \
#     --filter-name "ReadPosRankSum" \
#     --filter-expression "ReadPosRankSum < -8.0"

# ## Keep only SNPs that pass the above filters
# gatk SelectVariants \
#     -R $ref \
#     -V simchr_genotype_output_nobaseQrecal_filtered.vcf \
#     --select-type-to-include SNP \
#     -select 'vc.isNotFiltered()' \
#     -O filtered_simchr_nobaseQrecal_SNPs.vcf
