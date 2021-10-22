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
# Set Variables
# -----------------------------------------------------------------------------

#
## Load variables from settings file
source /home/aubkbk001/roh_param_project/init_script_vars.sh

## Set step name
STEP=04_downsample
PREV_STEP=03_read_align

## Create variable pointing to user's home directory for this step
HOME_DIR=/home/${USER}/${PROJECT}/${STEP}/

## Set the input directory, which is the output directory from the previous
## step
INPUT_DIR=/scratch/${USER}/${PROJECT}/${PREV_STEP}/output/sorted_bam

## Create a working directory on /scratch
mkdir -p /scratch/${USER}/${PROJECT}/${STEP}
chmod -R 700 /scratch/${USER}
WORK_DIR=/scratch/${USER}/${PROJECT}/${STEP}

## Create output directory in working directory
mkdir ${WORK_DIR}/output 
chmod -R 700 ${WORK_DIR} 
OUTPUT_DIR=${WORK_DIR}/output 

## cd into working scratch directory
cd ${WORK_DIR}


# -----------------------------------------------------------------------------
# do the downsample
# -----------------------------------------------------------------------------

# Create a mapping of downsample fractions to coverage levels. Used to display
# coverage levels in file and directory names
declare -A COVERX
COVERX[1.0]=50x
COVERX[0.6]=30x
COVERX[0.3]=15x
COVERX[0.2]=10x
COVERX[0.1]=05x

# -----------------------------------------------------------------------------
# Create text files listing sample individuals for each population size we
# want to test. Create directories for saving downsampled BAM files
# -----------------------------------------------------------------------------

for population in 100 50 30
do
    # Randomly choose n individuals from original population
    # Save list to file (pop_{population}_sample_list.txt)
	cd ${INPUT_DIR}

	ls i*.bam| sort -R | tail -${population} > ../sample_list.txt
	sed 's/_genome_sorted.bam//g' ../sample_list.txt > ${OUTPUT_DIR}/sample_id_list_pop_${population}.txt
    rm ../sample_list.txt

	for coverage in 1.0 0.6 0.3 0.2 0.1
	do
        mkdir ${OUTPUT_DIR}/sample_pop${population}_cvg_${COVERX[$coverage]}
        # echo ${COVERX[$coverage]}
	done
done


# -----------------------------------------------------------------------------
# Copy output files to user's home directory. If an output directory exists
# in the user's home directory, rename the existing output directory, then 
# copy the new output directory.
# -----------------------------------------------------------------------------

if [ -d "$HOME_DIR/output" ]; 
then
  TIMESTAMP=`date "+%Y%m%d-%H%M%S"` 
  mv ${HOME_DIR}/output ${HOME_DIR}/output_${TIMESTAMP}
fi

cp -r ${OUTPUT_DIR} ${HOME_DIR}/output

