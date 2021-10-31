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

## Load variables from settings file
source /home/aubkbk001/roh_param_project/init_script_vars.sh

## Set step name
STEP=04_downsample
PREV_STEP=03_read_align

## Create variable pointing to user's home directory for this project and step
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

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------
module load samtools/1.11

# -----------------------------------------------------------------------------
# Do the downsampling.
# -----------------------------------------------------------------------------

# Create subsets of various sizes of the initial sample population
# for population in 100 50 30
for population in ${popN[@]}
do
    # Randomly choose n individuals from original population. Save list of 
    # individuals to a text file
	cd ${INPUT_DIR}
	ls i*.bam| sort -R | tail -${population} > ../sample_list.txt
	sed 's/_genome_sorted.bam//g' ../sample_list.txt > ${OUTPUT_DIR}/sample_id_list_pop_${population}.txt
    rm ../sample_list.txt

    # Downsample each individual in the population to the specified levels of
    # coverage 
	for i in $(seq 0 $cvgCnt)
 	do
        # Create output directory for each population size and coverage level
        # combination.
        POP_CVG_OUTPUT_DIR=${OUTPUT_DIR}/sample_pop_${population}_cvg_${cvgX[i]}
        mkdir ${POP_CVG_OUTPUT_DIR}

        # Downsample each individual
        while read -a line
        do
            samtools view -s ${cvgP[i]} -@ 19 \
              -o ${POP_CVG_OUTPUT_DIR}/${line[0]}_pop_${population}_cvg_${cvgX[i]}.bam \
              ${INPUT_DIR}/${line[0]}_genome_sorted.bam
    
        done < ${OUTPUT_DIR}/sample_id_list_pop_${population}.txt
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

