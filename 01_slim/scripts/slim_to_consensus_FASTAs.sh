#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - ##### queue      |
#   |    - ## CPU + ## Gb   |
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

##  Set username
USER=aubaxh002

## Set project name
PROJ=01_slim

## Create a directory on /scratch
mkdir /scratch/${USER}_${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}_${PROJ}/

##  Copy input files to scratch
cp /home/$USER/roh_param_project/$PROJ/input/* /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/

## --------------------------------
## Load modules
module load anaconda
module load bcftools


## --------------------------------
## Run SLiM
# mkdir slim_output_files

# slim \
#   -d NUM_GEN=1000 \
#   -d MUTATION_RATE=100e-8 \
#   -d REC_RATE=250e-8 \
#   -d CHR_LENGTH=30000000 \
#   -d POP_SIZE=100 \
#   -d SAMPLE_SIZE=30 \
#   -d "OUT_PATH='/slim_output_files/'" \
#   /home/$USER/roh_param_project/$PROJ/scripts/chrom_evo_and_sample.slim


## --------------------------------
## Format SLiM output for read simulation
cd slim_output_files/

## Compress output VCF
# bgzip final_pop.vcf

## Index compressed VCF
tabix final_pop.vcf.gz

## Split output VCF into sample-specific VCF files
cd ../
mkdir sample_vcf_files
bcftools +split -O z -o sample_vcf_files/ ./slim_output_files/final_pop.vcf.gz


## Randomly select sample VCF files for conversation to FASTAs
## >> Would be more efficient to process 100 and downsample later, but this keeps
## >> each set of samples an independent draw from the final population
cd sample_vcf_files

for a in 15 30 50 100
	do
	ls i*.vcf.gz | sort -R | tail -100 > ../vcf_file_list_n.$a.txt
	sed 's/.vcf.gz//g' ../vcf_file_list_n.$a.txt > ../sample_id_list_n.$a.txt
	done


## Convert sample VCF files to two separate haplotype FASTAs
cd ../
mkdir sample_fasta_files
cd sample_fasta_files/

for a in 15 30 50 100
	do
	while read -a line
		do
		tabix ../sample_vcf_files/${line[0]}.vcf.gz
	
		bcftools consensus --haplotype 1 \
		--fasta-ref ../slim_output_files/ancestral.fasta \
		--output ${line[0]}_1_n.$a.fasta \
		../sample_vcf_files/${line[0]}.vcf.gz
	
		bcftools consensus --haplotype 2 \
		--fasta-ref ../slim_output_files/ancestral.fasta \
		--output ${line[0]}_2_n.$a.fasta \
		../sample_vcf_files/${line[0]}.vcf.gz
	
		done < ../sample_id_list_n.$a.txt
	done