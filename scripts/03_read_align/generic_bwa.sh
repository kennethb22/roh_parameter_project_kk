#!/bin/bash
#
#   +-----------------------+
#   |  USE:                 |
#   |    - LARGE queue      |
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

##  Set username
USER=aubaxh002

## Set project name
PROJ=psmc

## Create a directory on /scratch
mkdir /scratch/${USER}_${PROJ}/

## Set permissions for directory
chmod 700 /scratch/${USER}_${PROJ}/

##  Copy input files to scratch
cp /home/$USER/$PROJ/input/* /scratch/${USER}_${PROJ}/

## cd into working scratch directory
cd /scratch/${USER}_${PROJ}/


## --------------------------------
## Load modules
module load bwa/0.7.12
module load samtools/1.11


## --------------------------------
## Index reference genomes
# bwa index ./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa
# bwa index ./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa
	

## --------------------------------
## Align reads to reference genomes
## D genome - stephensi
while read -a line
	do
	bwa mem -t 20 -M \
	./sra_downloads/GCA_019054845.1_ASM1905484v1_genomic.fa \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	> ./align_files/${line[0]}_Dgenome.bam

## C genome	- stephensi 
	bwa mem -t 20 -M \
	./sra_downloads/GCA_001984765.1_C.can_genome_v1.0_genomic.fa \
	./cleaned_reads/trimmed_paired_${line[0]}_1.fastq.gz \
	./cleaned_reads/trimmed_paired_${line[0]}_2.fastq.gz \
	> ./align_files/${line[0]}_Cgenome.bam
	done < d_stephensi_sra_list.txt
	
## --------------------------------
## Merge and sort BAM files
cd ./align_files/

## D genome - stephensi
samtools sort -@ 19 sorted_d_stephensi_Dgenome.bam SRR14572526_Dgenome.bam 

## C genome	- stephensi
samtools sort -@ 19 sorted_d_stephensi_Cgenome.bam SRR14572526_Cgenome.bam
