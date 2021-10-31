#!/bin/bash
#
#  Replace the USER name in this script with your username and
#  call your project whatever you want
#
#  This script must be made executable like this
#    chmod +x my_script
#
#  Submit this script to the queue with a command like this
#    run_script my_script.sh

##  Set username
USER=aubaxh002

## Set project name
PROJ=04b_plink

## Create a directory on /scratch
# mkdir /scratch/${USER}_${PROJ}/

## Set permissions for directory
# chmod 700 /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/


## --------------------------------
## Load module 
module load anaconda/3-2020.11


## --------------------------------
## Convert VCF to PLINK format(s)
plink \
--vcf /scratch/aubaxh002_03b_noBQSR_genomicsdb/all_samples_nobaseQrecal_filtered_SNPs_only.vcf.gz \
--allow-extra-chr \
--out gatk_more_stringent


## --------------------------------
## Run PLINK to ID ROHs, varying a couple of parameters to check effects on output:
## try out a few different values for # of heterozygous sites allowed within a window,
## default = 1
for a in 1 2 3 4
do
## try out a few different values for max inverse density (kb/var),
## default = at least 1 var / 50 kb
	for b in 10 50 100 
	do
		plink --bim gatk_more_stringent.bim \
		--bed gatk_more_stringent.bed \
		--fam gatk_more_stringent.fam \
		--allow-extra-chr --no-pheno \
		--homozyg-window-het $a \
		--homozyg-window-missing 5 \
		--homozyg-window-snp 50 \
		--homozyg-density $b \
		--homozyg-gap 100 \
		--homozyg-window-threshold 0.05 \
		--homozyg-snp 100 \
		--homozyg-kb 100 \
		--out ./output/gatk_more_stringent.$a.$b.roh
	done
done


## --------------------------------
## Copy results back to project output directory (in home)



























