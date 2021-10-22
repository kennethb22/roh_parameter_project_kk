#!/bin/bash
#
#  +---------+
#  | BATCHED |
#  +---------+
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

## cd into directory
cd /scratch/${USER}_${PROJ}/

## set reference FASTA sequencd
ref=/scratch/aubaxh002_assem_indexing/hifiasm_kangaroo_rat_6cells.p_ctg.fasta

## --------------------------------
## Load modules
# module load samtools/1.11
module load picard/1.79
module load gatk/4.1.4.0


## --------------------------------
## Mark duplicate reads - GATK tools ignore them, no need to remove, just flag
## Will need to edit I option to match BAM files names used.
## Can use a file with one BAM file name per line (list_of_bam_files.txt) to loop
## over them.
while read -a line
	do
	java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/MarkDuplicates.jar \
	I=sorted_${line[0]}_small_genome_rgroups.bam \
	O=${line[0]}_small_genome_rgroups_dupmarked.bam \
	M=${line[0]}_rgroups_markdup_metrics.txt \
	MAX_RECORDS_IN_RAM=250000

	## Fix mate information
	java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/FixMateInformation.jar \
	I=${line[0]}_small_genome_rgroups_dupmarked.bam \
	O=${line[0]}_small_genome_rgroups_dupmarked_fixmate.bam

	## Index BAM files
	java -jar /opt/asn/apps/picard_1.79/picard_1.79/picard-tools-1.79/BuildBamIndex.jar \
	I=${line[0]}_small_genome_rgroups_dupmarked_fixmate.bam

	## Run HaplotypeCaller in GVCF mode
	gatk HaplotypeCaller \
		-R $ref \
		-I ${line[0]}_small_genome_rgroups_dupmarked_fixmate.bam \
		-stand-call-conf 20.0 \
		--emit-ref-confidence GVCF \
		-O ${line[0]}_nobaseQrecal.g.vcf
done < /home/aubaxh002/sample_lists/list_of_bam_files.txt

