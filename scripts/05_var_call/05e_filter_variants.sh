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
#cd
#  I ran each combination of population and coverage as separate jobs.
#
#  To do - modify script to create scripts to run for each population/coverage
#  combinationmore
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
declare -a cvgX=(50x 30x 15x 10x 05x)

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

        # ----------------------------------------------------------------------
        # Flag Variants that do not pass filters
        # ----------------------------------------------------------------------
        SCRIPT_FILE=05e_fltr_p${population}_c_${cvgX[i]}.sh
        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_gvcfs.vcf
        FILTERED_FILE=${OUT_FILE}
        GENOTYPED_FILE=sample_pop_${population}_cvg_${cvgX[i]}_genotyped_gvcfs.vcf

        # start_logging "gatk VariantFiltration- ${OUT_FILE}"

        # # Do VariantFiltration

        # gatk VariantFiltration \
        #     -R ${REF_GENOME_FILE} \
        #     -V ${OUTPUT_DIR}/${GENOTYPED_FILE} \
        #     -O ${OUTPUT_DIR}/${OUT_FILE} \
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

        echo '#!/bin/bash' >${SCRIPT_FILE}
        echo 'module load gatk/4.1.4.0' >>${SCRIPT_FILE}
        printf "\n\n" >>${SCRIPT_FILE}
        echo 'gatk VariantFiltration \' >>${SCRIPT_FILE}
        printf "%s %s \\" "-R" ${REF_GENOME_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s \\" "-V" ${OUTPUT_DIR} ${GENOTYPED_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s \\" "-O" ${OUTPUT_DIR} ${OUT_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-name" "QD" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-expression" "QD < 2.0" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-name" "FS" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-expression" "FS > 40.0" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-name" "SOR" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-expression" "SOR > 5.0" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-name" "MQ" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-expression" "MQ < 20.0" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-name" "MQRankSum" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-expression" "MQRankSum < -3.0 || MQRankSum > 3.0" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" \\' "--filter-name" "ReadPosRankSum" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf '%s \"%s\" ' "--filter-expression" "ReadPosRankSum < -8.0" >>${SCRIPT_FILE}
        printf "\n\n" >>${SCRIPT_FILE}

        # # ----------------------------------------------------------------------
        # # Keep only SNPs that pass the above filters
        # # ----------------------------------------------------------------------

        OUT_FILE=sample_pop_${population}_cvg_${cvgX[i]}_filtered_SNPs.vcf

        # start_logging "gatk SelectVariants - ${OUT_FILE}"

        # # Do Select Variants

        # gatk SelectVariants \
        #     -R ${REF_GENOME_FILE} \
        #     -V ${OUTPUT_DIR}/${FILTERED_FILE} \
        #     --select-type-to-include SNP \
        #     -select 'vc.isNotFiltered()' \
        #     -O ${OUTPUT_DIR}/${OUT_FILE}

        # stop_logging

        echo 'gatk SelectVariants \' >>${SCRIPT_FILE}
        printf "%s %s \\" "-R" ${REF_GENOME_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s \\" "-V" ${OUTPUT_DIR} ${FILTERED_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s \\" "--select-type-to-include SNP" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s \'%s\' \\" "-select" "vc.isNotFiltered()" >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}
        printf "%s %s/%s " "-O" ${OUTPUT_DIR} ${OUT_FILE} >>${SCRIPT_FILE}
        printf "\n" >>${SCRIPT_FILE}

        chmod +x ${SCRIPT_FILE}

        run_script ${SCRIPT_FILE}

    done

done

sleep 10
# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

# source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
