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
#  In user home directory, make a backup of the .asc_queue file and replace it
#  with .asc_queue_gatk_hc file before running haplotype calling. The .asc_queue
#  file sets the default parameters used by run_script.
#
#  This script submits one queue job for each each file on which gatk Haplotype
#  Caller is run, and in order for it to work it depends on the values in
#  .asc_queue being set to the values saved in .asc_queue_gatk_hc. Restore the
#  original .asc_queue file after running this script.
#
#  TO DO:  Split this script into two separate scripts: one for add read group
#          information and one for haplotype calling.

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=05_var_call
PREV_STEP=04_downsample
SCRIPT=05a_snp_calling

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
sleep 5

module load picard/1.79
# module load gatk/4.1.4.0

# -----------------------------------------------------------------------------
# Run Haplotype caller on all sample files
# -----------------------------------------------------------------------------

# Iterate over levels of coverage ------------------------------------------
declare -a cvgX=(30x 50x)
cvgCnt=${#cvgX[@]}
let cvgCnt-=1

for i in $(seq 0 $cvgCnt); do

    # Create output directory for each coverage level.

    CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_cvg_${cvgX[i]}
    mkdir ${CVG_OUTPUT_DIR}

    # Set input directory for this coverage

    CVG_INPUT_DIR=${INPUT_DIR}/sample_cvg_${cvgX[i]}

    # Process each individual sample file ----------------------------------

    while read -a line; do

        # # ------------------------------------------------------------------
        # # Add read group information
        # # ------------------------------------------------------------------

        # IN_FILE=${line[0]}_cvg_${cvgX[i]}.bam
        OUT_FILE=${line[0]}_cvg_${cvgX[i]}_rgroups.bam
        RGROUPS_FILE=${OUT_FILE}
        # start_logging "picard/AddOrReplaceReadGroups - ${OUT_FILE}"

        # ## Do picard read groups

        # java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar \
        #     I=${CVG_INPUT_DIR}/${IN_FILE} \
        #     O=${CVG_OUTPUT_DIR}/${OUT_FILE} SORT_ORDER=coordinate RGID=group1 RGLB=lib1 RGPL=illumina RGSM=${line[0]} RGPU=unit4

        # stop_logging

        # # ------------------------------------------------------------------
        # # Index BAM files
        # # ------------------------------------------------------------------

        # OUT_FILE=${line[0]}_cvg_${cvgX[i]}_rgroups.bai
        # start_logging "picard/BuildBamIndex - ${OUT_FILE}"

        # ## Run picard BuildBamIndex

        # java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/BuildBamIndex.jar \
        #     I=${CVG_OUTPUT_DIR}/${RGROUPS_FILE}

        # stop_logging

        # ------------------------------------------------------------------
        # Run HaplotypeCaller in GVCF mode
        # ------------------------------------------------------------------

        OUT_FILE=${line[0]}_cvg_${cvgX[i]}.g.vcf
        start_logging "gatk Haplotype Caller - ${OUT_FILE}"

        # Run gatk HaplotypeCaller

        # gatk HaplotypeCaller \
        #     -R ${REF_GENOME_FILE} \
        #     -I ${CVG_OUTPUT_DIR}/${RGROUPS_FILE} \
        #     -stand-call-conf 20.0 \
        #     --emit-ref-confidence GVCF \
        #     -O ${CVG_OUTPUT_DIR}/${OUT_FILE}

        SCRIPT_FILE=hc_${line[0]}_cvg_${cvgX[i]}.sh

        echo '#!/bin/bash' >${SCRIPT_FILE}
        echo 'module load gatk/4.1.4.0' >>${SCRIPT_FILE}
        echo 'gatk HaplotypeCaller \' >>${SCRIPT_FILE}
        printf "%s %s \\" "-R" ${REF_GENOME_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s \\" "-I" ${CVG_OUTPUT_DIR}/${RGROUPS_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        echo '-stand-call-conf 20.0 \' >>${SCRIPT_FILE}
        echo '--emit-ref-confidence GVCF \' >>${SCRIPT_FILE}
        printf "%s %s/%s \\" "-O" ${CVG_OUTPUT_DIR}/${OUT_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}

        chmod +x ${SCRIPT_FILE}

        run_script ${SCRIPT_FILE}

        sleep 10
        # stop_logging

    done <${SAMPLE_ID_LIST}
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
