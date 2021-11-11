#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - ?LARGE queue     |
#   |    - 20 CPU + 20 Gb   |
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

# -----------------------------------------------------------------------------
# Set variables for this step
# -----------------------------------------------------------------------------

STEP=04_downsample
PREV_STEP=03_read_align
SCRIPT=04_run_downsample.sh

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load samtools/1.11

# -----------------------------------------------------------------------------
# Do the downsampling.
# -----------------------------------------------------------------------------

# Create subsets of various sizes of the initial sample population
# for population in 100 50 30

for population in ${popN[@]}; do

    # Randomly choose n individuals from original population. Save list of
    # individuals to a text file

    cd ${INPUT_DIR}

    ls i*_sorted.bam | sort -R | tail -${population} >../sample_list.txt
    sed 's/_genome_sorted.bam//g' ../sample_list.txt >${OUTPUT_DIR}/sample_id_list_pop_${population}.txt
    rm ../sample_list.txt

    # Downsample each individual in the population to the specified levels of
    # coverage

    for i in $(seq 0 $cvgCnt); do

        # Create output directory for each population size and coverage level
        # combination.

        POP_CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}
        mkdir ${POP_CVG_OUTPUT_DIR}

        # Downsample each individual

        while read -a line; do

            OUT_FILE=${line[0]}_pop_${population}_cvg_${cvgX[i]}.bam
            start_logging "samtools view - ${OUT_FILE}"

            samtools view -s ${cvgP[i]} -@ 19 \
                -o ${POP_CVG_OUTPUT_DIR}/${OUT_FILE} \
                ${INPUT_DIR}/${line[0]}_genome_sorted.bam

            stop_logging

        done <${OUTPUT_DIR}/sample_id_list_pop_${population}.txt
    done
done

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
