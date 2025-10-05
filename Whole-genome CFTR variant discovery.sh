#!/bin/bash
# Whole-genome CFTR variant discovery pipeline (child + father)


# Step 1: QC
fastqc child_1.fastq.gz child_2.fastq.gz father_1.fastq.gz father_2.fastq.gz -o qc_reports/


# Step 2: Alignment
bwa mem -R "@RG\tID:child\tSM:child\tPL:ILLUMINA" Homo_sapiens_assembly38.fasta \
child_1.fastq.gz child_2.fastq.gz | samtools sort -o child.sorted.bam


bwa mem -R "@RG\tID:father\tSM:father\tPL:ILLUMINA" Homo_sapiens_assembly38.fasta \
father_1.fastq.gz father_2.fastq.gz | samtools sort -o father.sorted.bam


# Step 3: Mark Duplicates
gatk MarkDuplicates -I child.sorted.bam -O child.markdup.bam -M child.metrics.txt
gatk MarkDuplicates -I father.sorted.bam -O father.markdup.bam -M father.metrics.txt


# Step 4: BQSR
gatk BaseRecalibrator -I child.markdup.bam -R Homo_sapiens_assembly38.fasta \
--known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
--known-sites Homo_sapiens_assembly38.known_indels.vcf.gz -O child.recal.table


gatk ApplyBQSR -R Homo_sapiens_assembly38.fasta -I child.markdup.bam --bqsr-recal-file child.recal.table -O child.final.bam


# Repeat for father

gatk BaseRecalibrator -I father.markdup.bam -R Homo_sapiens_assembly38.fasta \
--known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
--known-sites Homo_sapiens_assembly38.known_indels.vcf.gz -O father.recal.table

gatk ApplyBQSR -R Homo_sapiens_assembly38.fasta -I father.markdup.bam --bqsr-recal-file father.recal.table -O father.final.bam



# Step 5: Variant Calling
gatk HaplotypeCaller -R Homo_sapiens_assembly38.fasta -I child.final.bam -O child.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R Homo_sapiens_assembly38.fasta -I father.final.bam -O father.g.vcf.gz -ERC GVCF


# Step 6: Joint Genotyping
gatk CombineGVCFs -R Homo_sapiens_assembly38.fasta -V child.g.vcf.gz -V father.g.vcf.gz -O family.g.vcf.gz
gatk GenotypeGVCFs -R Homo_sapiens_assembly38.fasta -V family.g.vcf.gz -O family.raw.vcf.gz


# Step 7: Filtering
gatk VariantFiltration -R Homo_sapiens_assembly38.fasta -V family.raw.vcf.gz \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "FAIL" -O family.filtered.vcf.gz


# Step 8: Extract CFTR locus
tabix -h family.filtered.vcf.gz 7:117120016-117308718 > CFTR_family.vcf


# Step 9: Annotate with VEP
vep -i CFTR_family.vcf -o CFTR_family_annotated.vcf --cache --assembly GRCh38 --vcf --everything


# Step 10: Phase with WhatsHap
whatshap phase --reference Homo_sapiens_assembly38.fasta --ped family.ped CFTR_family.vcf child.final.bam -o CFTR_family.phased.vcf
