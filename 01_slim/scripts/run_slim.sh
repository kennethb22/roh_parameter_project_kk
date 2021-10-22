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


## Load variables from settings file
source /home/aubkbk001/roh_param_project/init_script_vars.sh

## Set step name
STEP=01_slim

## Create variable pointing to user's home directory for this step
HOME_DIR=/home/${USER}/${PROJECT}/${STEP}/

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


## --------------------------------
## Run SLiM
module load anaconda/3-2020.11

# final run settings (w/ 10e3 generations):
for m in 5e-07
	do
	for r in 1e-8
		do
		for p in  500 
			do
			mkdir ${OUTPUT_DIR}/slim_output_files_m${m}_r${r}_p${p}
				# models chromosome with coding and non-coding regions
				slim \
					-d POP_SIZE=$p \
					-d MUTATION_RATE=$m \
					-d RECOMB_RATE=$r \
					-d "OUT_PATH='${OUTPUT_DIR}/slim_output_files_m${m}_r${r}_p${p}/'" \
					/home/$USER/$PROJECT/$STEP/scripts/chrom_w_struct_and_evo.slim
			done
		done
	done

mail -s 'SLiM run finished - Starting BCFtools' ${EMAIL} <<< 'SLiM run finished - Starting BCFtools'

# ## --------------------------------
# ## Format SLiM output for read simulation
# module purge
# module load bcftools/1.13

# ## final run settings (w/ 10e3 generations):
# for m in 5e-07
# 	do
# 	for r in 1e-8
# 		do
# 		for p in  500 
# 			do
# 			cd /scratch/${USER}_${PROJ}/slim_output_files_m${m}_r${r}_p${p}/

# 			## Compress output VCF
# 			bgzip final_pop.vcf

# 			## Index compressed VCF
# 			tabix final_pop.vcf.gz

# 			## Split output VCF into sample-specific VCF files
# 			cd ../
# 			mkdir sample_vcf_files_m${m}_r${r}_p${p}
# 			bcftools +split -O z -o sample_vcf_files_m${m}_r${r}_p${p}/ ./slim_output_files_m${m}_r${r}_p${p}/final_pop.vcf.gz


# 			## Randomly select sample VCF files for conversation to FASTAs
# 			## >> Selecting 100, can be downsampled later
# 			cd sample_vcf_files_m${m}_r${r}_p${p}/

# 			ls i*.vcf.gz | sort -R | tail -100 > ../vcf_file_list_m${m}_r${r}_p${p}.txt
# 			sed 's/.vcf.gz//g' ../vcf_file_list_m${m}_r${r}_p${p}.txt > ../sample_id_list_m${m}_r${r}_p${p}.txt


# 			## Convert sample VCF files to two separate haplotype FASTAs per individual
# 			cd ../
# 			mkdir sample_fasta_files_m${m}_r${r}_p${p}
# 			cd sample_fasta_files_m${m}_r${r}_p${p}/

# 			while read -a line
# 				do
# 				tabix ../sample_vcf_files_m${m}_r${r}_p${p}/${line[0]}.vcf.gz
	
# 				bcftools norm --check-ref s \
# 				--fasta-ref ../slim_output_files_m${m}_r${r}_p${p}/ancestral.fasta \
# 				--multiallelics - \
# 				--do-not-normalize \
# 				--output ../sample_vcf_files_m${m}_r${r}_p${p}/norm_${line[0]}.vcf \
# 				../sample_vcf_files_m${m}_r${r}_p${p}/${line[0]}.vcf.gz
	
# 				bgzip ../sample_vcf_files_m${m}_r${r}_p${p}/norm_${line[0]}.vcf
# 				tabix ../sample_vcf_files_m${m}_r${r}_p${p}/norm_${line[0]}.vcf.gz

# 				bcftools +fixref ../sample_vcf_files_m${m}_r${r}_p${p}/norm_${line[0]}.vcf.gz -- \
# 				--fasta-ref ../slim_output_files_m${m}_r${r}_p${p}/ancestral.fasta
		
# 				bcftools consensus --haplotype 1 \
# 				--fasta-ref ../slim_output_files_m${m}_r${r}_p${p}/ancestral.fasta \
# 				--output ${line[0]}_1.fasta \
# 				../sample_vcf_files_m${m}_r${r}_p${p}/norm_${line[0]}.vcf.gz

# 				bcftools consensus --haplotype 2 \
# 				--fasta-ref ../slim_output_files_m${m}_r${r}_p${p}/ancestral.fasta \
# 				--output ${line[0]}_2.fasta \
# 				../sample_vcf_files_m${m}_r${r}_p${p}/norm_${line[0]}.vcf.gz

# 				done < ../sample_id_list_m${m}_r${r}_p${p}.txt
# 			done
# 		done
# 	done


# mail -s 'SLiM run finished - submit Rviz' ${EMAIL} <<< 'SLiM run finished'

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
