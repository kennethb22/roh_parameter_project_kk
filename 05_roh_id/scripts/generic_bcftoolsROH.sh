#!/bin/bash
#
#  +-----------------+
#  | REQUIRES 20 CPU |
#  +-----------------+
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

##  Set username
USER=aubaxh002

## Set project name
PROJ=04a_bcftoolsROH

## Create a directory on /scratch
# mkdir /scratch/${USER}_${PROJ}/

## Set permissions for directory
# chmod 700 /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/


## --------------------------------
## Load modules
module load bcftools/1.10.2
module load samtools/1.11


## --------------------------------
## Generate allele frequency files for use with bcftools roh & index with tabix
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' \
/scratch/aubaxh002_03b_noBQSR_genomicsdb/all_samples_nobaseQrecal_filtered_SNPs_only.vcf.gz \
| bgzip -c > freqs_morestringent_filteredSNPsonly.tab.gz

tabix -s1 -b2 -e2 freqs_morestringent_filteredSNPsonly.tab.gz

## --------------------------------
## Try ROH calling with AF estimation using GT, compare with PL estimation output

## SNPs filtered using more stringent criteria
bcftools roh \
--GTs-only 30 \
--threads 20 \
-o all_samples_GTonly_AFest_roh_morestringent.txt \
--AF-file freqs_morestringent_filteredSNPsonly.tab.gz \
/scratch/aubaxh002_03b_noBQSR_genomicsdb/all_samples_nobaseQrecal_filtered_SNPs_only.vcf.gz

bcftools roh \
--threads 20 \
-o all_samples_GTPL_AFest_roh_morestringent.txt \
--AF-file freqs_morestringent_filteredSNPsonly.tab.gz \
/scratch/aubaxh002_03b_noBQSR_genomicsdb/all_samples_nobaseQrecal_filtered_SNPs_only.vcf.gz


## --------------------------------
## Extract information on ROHs (i.e., exclude information on individual sites)
grep "^RG" all_samples_GTonly_AFest_roh_morestringent.txt >\
all_samples_GTonly_AFest_roh_morestringent_RG_ONLY.txt

grep "^RG" all_samples_GTPL_AFest_roh_morestringent.txt >\
all_samples_GTPL_AFest_roh_morestringent_RG_ONLY.txt

