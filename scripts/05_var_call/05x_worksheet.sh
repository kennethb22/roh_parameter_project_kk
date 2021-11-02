#!/bin/bash
#
#  +----------------------+
#  | REQUEST 4 CPU + 60gb |
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

##  Set username
USER=aubaxh002

## Set project name
PROJ=04_var_call

## Create a directory on /scratch
mkdir /scratch/${USER}_${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/

## --------------------------------
## Load modules
module load gatk/4.1.4.0
module load picard/1.79
module load vcftools/0.1.14

## --------------------------------
## Genotype GVCFs across all samples simultaneously
## reference FASTA sequence file
ref=/scratch/aubaxh002_assem_indexing/hifiasm_kangaroo_rat_6cells.p_ctg.fasta

## -V [sample VCF file]; unfortunately, no easy way to just provide a file with names;
## will need one -V line per sample

OUT_DIR=/scratch/aubkbk001/roh_param_project/data/05_var_call/output/sample_pop_03_cvg_05x

gatk GenomicsDBImport \
    -V ${OUT_DIR}/i31_pop_03_cvg_05x.g.vcf \
    -V ${OUT_DIR}/i322_pop_03_cvg_05x.g.vcf \
    -V ${OUT_DIR}/i428_pop_03_cvg_05x.g.vcf \
    -L 1 \
    --genomicsdb-workspace-path simchr_database_gatk_genomics

gatk CombineGVCFs \
    -R /scratch/aubkbk001/roh_param_project/data/01_slim/output/slim_output_files_m5e-07_r1e-8_p500/ancestral.fasta \
    -V ${OUT_DIR}/i31_pop_03_cvg_05x.g.vcf \
    -V ${OUT_DIR}/i322_pop_03_cvg_05x.g.vcf \
    -V ${OUT_DIR}/i428_pop_03_cvg_05x.g.vcf \
    -O cohort.g.vcf.gz

gatk ValidateVariants \
    -V ${OUT_DIR}/i31_pop_03_cvg_05x.g.vcf \
    -R /scratch/aubkbk001/roh_param_project/data/01_slim/output/slim_output_files_m5e-07_r1e-8_p500/ancestral.fasta

gatk GenotypeGVCFs \
    -R $ref \
    -V gendb://simchr_database_gatk_genomics \
    -O simchr_genotype_output_nobaseQrecal.vcf

## Flag variants that do not pass filters
gatk VariantFiltration \
    -R $ref \
    -V simchr_genotype_output_nobaseQrecal.vcf \
    -O simchr_genotype_output_nobaseQrecal_filtered.vcf \
    --filter-name "QD" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS" \
    --filter-expression "FS > 40.0" \
    --filter-name "SOR" \
    --filter-expression "SOR > 5.0" \
    --filter-name "MQ" \
    --filter-expression "MQ < 20.0" \
    --filter-name "MQRankSum" \
    --filter-expression " MQRankSum < -3.0 || MQRankSum > 3.0" \
    --filter-name "ReadPosRankSum" \
    --filter-expression "ReadPosRankSum < -8.0"

## Keep only SNPs that pass the above filters
gatk SelectVariants \
    -R $ref \
    -V simchr_genotype_output_nobaseQrecal_filtered.vcf \
    --select-type-to-include SNP \
    -select 'vc.isNotFiltered()' \
    -O filtered_simchr_nobaseQrecal_SNPs.vcf
