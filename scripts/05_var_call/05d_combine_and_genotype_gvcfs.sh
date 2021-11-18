#!/bin/bash
#
#  +----------------------+
#  | large queue:  REQUEST 2 CPU + 60 gb |
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
#  I ran each combination of population and coverage as separate jobs.
#
#  To do - modify script to create scripts to run for each population/coverage
#  combination
#
#  NOTES:
#
#  2021-11-14 - 05d_genotype_and_filter_p100_10x completed but did't
#    produce filtered.vcf output file. Checked job output log, execution just
#    stopped during the CombineVCFS step. I think that happened because I had
#    exceeded my disk quota. Scripts were making a new backup of the output
#    directory from scratch for every copy of 05d being run. Filled up quota
#    quickly. Disabling backup code in the script. Will manually backup of
#    output directory from scratch to home after all 05d scripts finish.
#

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
# Load Modules
# -----------------------------------------------------------------------------

module load gatk/4.1.4.0

# -----------------------------------------------------------------------------
# Genotype GVCFs across all samples simultaneously
# -----------------------------------------------------------------------------

## Create arrays of downsample levels to be used.
## Coverage level in NNx for display in output file names
# declare -a cvgX=(50x 30x 15x 10x 05x)
declare -a cvgX=(05x)
## Coverage level fraction to supply to samtools
declare -a cvgP=(1.0 0.6 0.3 0.2 0.1)

## Get length of the coverage level arrays. Subtract 1 because arrays are zero
## based, and we'll iterate over the arrays from 0 to cvgCnt
cvgCnt=${#cvgX[@]}
let cvgCnt-=1

## Create array of population sizes we want to test
# declare -a popN=(100 50 30)
declare -a popN=(100 50 30)

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

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_combined_gvcfs.vcf
        COMBINED_GVCFS_FILE=${OUT_FILE}
        SCRIPT_FILE=05d_gtyp_p${population}_c_${cvgX[i]}.sh

        echo '#!/bin/bash' >${SCRIPT_FILE}
        echo 'module load gatk/4.1.4.0' >>${SCRIPT_FILE}
        printf "\n\n" >>${SCRIPT_FILE}
        echo 'gatk CombineGVCFs \' >>${SCRIPT_FILE}
        printf "%s %s \\" "-R" ${REF_GENOME_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s \\" "-V" ${MAP_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s" "-O" ${OUTPUT_DIR} ${OUT_FILE} >>${SCRIPT_FILE}
        printf "\n\n" >>${SCRIPT_FILE}

        # ----------------------------------------------------------------------
        # gatk GenotypeGVCFs
        # ----------------------------------------------------------------------

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_genotyped_gvcfs.vcf
        GENOTYPED_FILE=${OUT_FILE}

        echo 'gatk GenotypeGVCFs \' >>${SCRIPT_FILE}
        printf "%s %s \\" "-R" ${REF_GENOME_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s \\" "-V" ${OUTPUT_DIR} ${COMBINED_GVCFS_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s" "-O" ${OUTPUT_DIR} ${OUT_FILE} >>${SCRIPT_FILE}
        printf "\n\n" >>${SCRIPT_FILE}

        chmod +x ${SCRIPT_FILE}

        run_script ${SCRIPT_FILE}
    done

done

# sleep 10
# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
