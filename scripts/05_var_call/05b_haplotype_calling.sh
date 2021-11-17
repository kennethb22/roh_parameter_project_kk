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
SCRIPT=$(echo "$(echo "$0" | sed -e "s/^\.\///")" | sed -e "s/\.sh//")

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load gatk/4.1.4.0

# -----------------------------------------------------------------------------
# Run Haplotype caller on all sample files
# -----------------------------------------------------------------------------

# Iterate over levels of coverage ------------------------------------------
declare -a cvgX=(30x)
cvgCnt=${#cvgX[@]}
let cvgCnt-=1

for i in $(seq 0 $cvgCnt); do

    # Create output directory for each coverage level.

    CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_cvg_${cvgX[i]}
    mkdir ${CVG_OUTPUT_DIR}

    # Process each individual sample file ----------------------------------

    while read -a line; do

        # ------------------------------------------------------------------
        # Run HaplotypeCaller in GVCF mode
        # ------------------------------------------------------------------

        OUT_FILE=${line[0]}_cvg_${cvgX[i]}.g.vcf
        RGROUPS_FILE=${line[0]}_cvg_${cvgX[i]}_rgroups.bam

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
        stop_logging

    done <${SAMPLE_ID_LIST}
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
