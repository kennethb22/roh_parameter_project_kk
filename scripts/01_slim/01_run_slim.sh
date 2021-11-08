#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - MEDIUM queue     |
#   |    - 2 CPU + 6 Gb     |
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

STEP=01_slim
PREV_STEP=05_var_call
SCRIPT=01_run_slim

# -----------------------------------------------------------------------------
# Load variables and functions from settings file
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/init_script_vars.sh

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module load anaconda/3-2020.11

# -----------------------------------------------------------------------------
# Run SLiM
# -----------------------------------------------------------------------------

SLIM_PARAM_FILE=/home/${USER}/${PROJECT}/scripts/${STEP}/chrom_w_struct_and_evo.slim

for m in 5e-07; do
    for r in 1e-8; do
        for p in 500; do
            SLIM_OUT_DIR=slim_m${m}_r${r}_p${p}
            start_logging "Run SLiM - ${SLIM_OUT_DIR}"
            mkdir ${OUTPUT_DIR}/${SLIM_OUT_DIR}

            # models chromosome with coding and non-coding regions
            slim \
                -d POP_SIZE=$p \
                -d MUTATION_RATE=$m \
                -d RECOMB_RATE=$r \
                -d "OUT_PATH='${OUTPUT_DIR}/${SLIM_OUT_DIR}/'" \
                ${SLIM_PARAM_FILE}
        done
    done
done

stop_logging

# mail -s 'SLiM run finished - Starting BCFtools' ${EMAIL} <<<'SLiM run finished - Starting BCFtools'

# -----------------------------------------------------------------------------
# Format SLiM output for read simulation
# -----------------------------------------------------------------------------

module purge
module load bcftools/1.13

start_logging "Format SLiM Output for read simulation - ${SLIM_OUT_DIR}"

for m in 5e-07; do
    for r in 1e-8; do
        for p in 500; do

            FASTA_OUT_DIR=sample_fasta_files_m${m}_r${r}_p${p}
            VCF_OUT_DIR=sample_vcf_files_m${m}_r${r}_p${p}

            cd ${OUTPUT_DIR}/${SLIM_OUT_DIR}

            ## Compress output VCF
            bgzip final_pop.vcf

            ## Index compressed VCF
            tabix final_pop.vcf.gz

            ## Split output VCF into sample-specific VCF files
            cd ../
            mkdir ${VCF_OUT_DIR}
            bcftools +split -O z -o ${VCF_OUT_DIR}/ ./${SLIM_OUT_DIR}/final_pop.vcf.gz

            ## Randomly select sample VCF files for conversation to FASTAs
            ## >> Selecting 100, can be downsampled later
            cd ${VCF_OUT_DIR}/

            ls i*.vcf.gz | sort -R | tail -100 >../vcf_file_list_m${m}_r${r}_p${p}.txt
            sed 's/.vcf.gz//g' ../vcf_file_list_m${m}_r${r}_p${p}.txt >../sample_id_list_m${m}_r${r}_p${p}.txt

            ## Convert sample VCF files to two separate haplotype FASTAs per individual
            cd ../
            mkdir ${FASTA_OUT_DIR}
            cd ${FASTA_OUT_DIR}/

            while read -a line; do
                tabix ../${VCF_OUT_DIR}/${line[0]}.vcf.gz

                bcftools norm --check-ref s \
                    --fasta-ref ../${SLIM_OUT_DIR}/ancestral.fasta \
                    --multiallelics - \
                    --do-not-normalize \
                    --output ../${VCF_OUT_DIR}/norm_${line[0]}.vcf \
                    ../${VCF_OUT_DIR}/${line[0]}.vcf.gz

                bgzip ../${VCF_OUT_DIR}/norm_${line[0]}.vcf
                tabix ../${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

                bcftools +fixref ../${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz -- \
                    --fasta-ref ../${SLIM_OUT_DIR}/ancestral.fasta

                bcftools consensus --haplotype 1 \
                    --fasta-ref ../${SLIM_OUT_DIR}/ancestral.fasta \
                    --output ${line[0]}_1.fasta \
                    ../${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

                bcftools consensus --haplotype 2 \
                    --fasta-ref ../${SLIM_OUT_DIR}/ancestral.fasta \
                    --output ${line[0]}_2.fasta \
                    ../${VCF_OUT_DIR}/norm_${line[0]}.vcf.gz

            done <../sample_id_list_m${m}_r${r}_p${p}.txt
        done
    done
done

stop_logging

# mail -s 'SLiM run finished - submit Rviz' ${EMAIL} <<<'SLiM run finished'

# -----------------------------------------------------------------------------
# Copy output files to user's home directory.
# -----------------------------------------------------------------------------

source /home/aubkbk001/roh_param_project/scripts/99_includes/backup_output.sh
