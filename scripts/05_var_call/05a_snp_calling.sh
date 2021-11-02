#!/bin/bash
#
#  +---------+
#  | BATCHED |
#  +---------+
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
SCRIPT=05a_snp_calling

## Create variable pointing to user's home directory for this project and step
HOME_STEP_DIR=/home/${USER}/${PROJECT}/data/${STEP}/

## Set the input directory, which is the output directory from the previous
## step
INPUT_DIR=/scratch/${USER}/${PROJECT}/data/${PREV_STEP}/output

## Create a working directory on /scratch
mkdir -p /scratch/${USER}/${PROJECT}/data/${STEP}
chmod -R 700 /scratch/${USER}
WORK_DIR=/scratch/${USER}/${PROJECT}/data/${STEP}

## Create output directory in working directory
mkdir ${WORK_DIR}/output
chmod -R 700 ${WORK_DIR}
OUTPUT_DIR=${WORK_DIR}/output

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load picard/1.79
module load gatk/4.1.4.0

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
# Run Haplotype caller on all sample files
# -----------------------------------------------------------------------------

# Iterate over each combination of population and coverage level. Sets of
# population and coverage are set in init_script_vars.sh

# Iterate over populations -----------------------------------------------------

for population in ${popN[@]}; do

    # Iterate over levels of coverage ------------------------------------------

    for i in $(seq 0 $cvgCnt); do

        # Create output directory for each population size and coverage level
        # combination.

        POP_CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}
        mkdir ${POP_CVG_OUTPUT_DIR}

        # Set input directory for this population and coverage

        POP_CVG_INPUT_DIR=${INPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}

        # Set name of map file for this population and coverage. This file
        # is used in 05b_genomicsdb.sh

        MAP_FILE=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}_map.list

        # Delete map file if it exists

        if [ -f "$MAP_FILE" ]; then
            rm ${MAP_FILE}
        fi

        # Process each individual sample file ----------------------------------

        while read -a line; do

            # Add read group information --------------------------------------

            # ## Set action and start time for log

            # OUT_FILE=${line[0]}_pop_${population}_cvg_${cvgX[i]}_rgroups.bam
            # ACTION="picard/AddOrReplaceReadGroups - ${OUT_FILE}"
            # START_TIME_HHMM=$(date +" %T")
            # START_TIME=$(date +%s)

            # ## Do picard read groups

            # java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar \
            #     I=${POP_CVG_INPUT_DIR}/${line[0]}_pop_${population}_cvg_${cvgX[i]}.bam \
            #     O=$${POP_CVG_OUTPUT_DIR}/{OUT_FILE} SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=${line[0]} RGPU=unit4

            # ## Write line to log file for this action

            # END_TIME_HHMM=$(date +"%T")
            # END_TIME=$(date +%s)
            # ELAPSED=$(expr $END_TIME - $START_TIME)
            # DURATION=$(date -u -d @${ELAPSED} +"%T")
            # printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}

            # ## Index BAM files -------------------------------------------------

            # ## Set action and start time for log

            # OUT_FILE=${line[0]}_pop_${population}_cvg_${cvgX[i]}_rgroups.bai
            # ACTION="picard/BuildBamIndex - ${OUT_FILE}"
            # START_TIME_HHMM=$(date +" %T")
            # START_TIME=$(date +%s)

            # ## Run picard BuildBamIndex

            # java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/BuildBamIndex.jar \
            #     I=${POP_CVG_OUTPUT_DIR}/${line[0]}_pop_${population}_cvg_${cvgX[i]}_rgroups.bam

            # ## Write line to log file for this action

            # END_TIME_HHMM=$(date +"%T")
            # END_TIME=$(date +%s)
            # ELAPSED=$(expr $END_TIME - $START_TIME)
            # DURATION=$(date -u -d @${ELAPSED} +"%T")
            # printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}

            # ## Run HaplotypeCaller in GVCF mode --------------------------------

            # ## Set action and start time for log

            # OUT_FILE=${line[0]}_pop_${population}_cvg_${cvgX[i]}.g.vcf
            # ACTION="gatk Haplotype Caller - ${OUT_FILE}"
            # START_TIME_HHMM=$(date +" %T")
            # START_TIME=$(date +%s)

            # ## Run gatk HaplotypeCaller

            # gatk HaplotypeCaller \
            #     -R ${REF_GENOME_FILE} \
            #     -I ${POP_CVG_OUTPUT_DIR}/${line[0]}_pop_${population}_cvg_${cvgX[i]}_rgroups.bam \
            #     -stand-call-conf 20.0 \
            #     --emit-ref-confidence GVCF \
            #     -O ${POP_CVG_OUTPUT_DIR}/${OUT_FILE}

            # ## Write line to log file for this action

            # END_TIME_HHMM=$(date +"%T")
            # END_TIME=$(date +%s)
            # ELAPSED=$(expr $END_TIME - $START_TIME)
            # DURATION=$(date -u -d @${ELAPSED} +"%T")
            # printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}

            # Write line for sample to map file.

            printf "%s\t%s\n" \
                ${POP_CVG_OUTPUT_DIR}/${line[0]}_pop_${population}_cvg_${cvgX[i]}.g.vcf \
                >>$MAP_FILE

        done <${INPUT_DIR}/sample_id_list_pop_${population}.txt
    done
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
