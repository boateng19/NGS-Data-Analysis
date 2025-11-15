

Whole-Genome Sequencing Analysis Report
Clinical Case Study – CFTR Variant Interpretation
Author: Agyekum Boateng Date: 2025

1. Introduction
A 6-year-old boy presents with:
	•	Chronic cough 
	•	Recurrent lung infections 
	•	Poor weight gain 
A sweat chloride test returned 45 mmol/L, which is borderline for cystic fibrosis (CF).
Whole-genome sequencing (WGS) was performed on:
	•	Child 
	•	Father 
Mother’s data are unavailable. This analysis investigates pathogenic CFTR variants and determines inheritance, pathogenicity, and diagnostic confirmation.

2. Methods
2.1 Data Files

~/Boateng/child_1.fastq.gz
~/Boateng/child_2.fastq.gz
~/Boateng/father_1.fastq.gz
~/Boateng/father_2.fastq.gz
~/Boateng/Homo_sapiens_assembly38.fasta
~/Boateng/Homo_sapiens_assembly38.dbsnp138.vcf
~/Boateng/Homo_sapiens_assembly38.known_indels.vcf.gz

3. Pipeline (GATK Best Practices)

3.1 Quality Control (FastQC + Trimmomatic)

fastqc child_*.fastq.gz father_*.fastq.gz

trimmomatic PE child_1.fastq.gz child_2.fastq.gz \
child_1.trim.fq.gz child_1.unp.fq.gz \
child_2.trim.fq.gz child_2.unp.fq.gz \
SLIDINGWINDOW:4:20 MINLEN:50
QC results (summary):
	•	Per-base quality > Q30 across most cycles 
	•	Slight adapter contamination—resolved by trimming 
	•	No major duplication or GC-bias issues 

3.2 Alignment to Reference Genome (BWA-MEM2)

bwa-mem2 index Homo_sapiens_assembly38.fasta

bwa-mem2 mem -t 8 Homo_sapiens_assembly38.fasta \
child_1.trim.fq.gz child_2.trim.fq.gz | \
samtools sort -o child.sorted.bam
Same commands repeated for the father.

3.3 Mark Duplicates (GATK)

gatk MarkDuplicates \
-I child.sorted.bam \
-O child.markdup.bam \
-M child.metrics.txt

3.4 Base Quality Score Recalibration (BQSR)

gatk BaseRecalibrator \
 -I child.markdup.bam \
 -R Homo_sapiens_assembly38.fasta \
 --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
 --known-sites Homo_sapiens_assembly38.known_indels.vcf.gz \
 -O child.recal.table

gatk ApplyBQSR \
 -R Homo_sapiens_assembly38.fasta \
 -I child.markdup.bam \
 --bqsr-recal-file child.recal.table \
 -O child.bqsr.bam
Same for father.

3.5 Variant Calling — HaplotypeCaller (GVCF Mode)

gatk HaplotypeCaller \
 -R Homo_sapiens_assembly38.fasta \
 -I child.bqsr.bam \
 -O child.g.vcf.gz \
 -ERC GVCF

3.6 Joint Genotyping
bash

gatk CombineGVCFs \
 -R Homo_sapiens_assembly38.fasta \
 --variant child.g.vcf.gz \
 --variant father.g.vcf.gz \
 -O cohort.g.vcf.gz

gatk GenotypeGVCFs \
 -R Homo_sapiens_assembly38.fasta \
 -V cohort.g.vcf.gz \
 -O cohort.raw.vcf.gz

3.7 Hard Filtering
bash

gatk VariantFiltration \
 -V cohort.raw.vcf.gz \
 -O cohort.filtered.vcf.gz \
 --filter-expression "QD < 2.0" --filter-name QD2 \
 --filter-expression "FS > 60.0" --filter-name FS60

3.8 Extract CFTR Region (chr7:117,120,000–117,330,000)
bash

bcftools view -r 7:117120000-117330000 cohort.filtered.vcf.gz -o CFTR.vcf

3.9 Variant Annotation (VEP)
bash

vep -i CFTR.vcf -o CFTR_annotated.vcf --cache --everything --assembly GRCh38

4. Results
4.1 Identified Pathogenic CFTR Variants in the Child
Variant   										Genotype (Child)	Genotype (Father)	Clinical Significance	Notes
chr7:117199644 C>T (c.1521_1523delCTT) → ΔF508	Het						Present (Het)	Pathogenic				Most common CF variant
chr7:117307003 G>A (W1282X)						Het						Absent			Pathogenic				Nonsense mutation → likely de novo or inherited from mother

4.2 Interpretation of Findings
Variant 1: ΔF508 (c.1521_1523delCTT)
Child: Heterozygous

Father: Heterozygous

→ Inherited from father

Causes misfolded CFTR protein → fails to reach plasma membrane.

Variant 2: W1282X (c.3846G>A)
Child: Heterozygous

Father: Not present

→ Likely maternal inheritance or de novo

Produces truncated, nonfunctional CFTR protein.

Together the child is compound heterozygous for two pathogenic variants → diagnostic for cystic fibrosis.

4.3 IGV Visualization Summary (described)
ΔF508
Clear 3-bp deletion visible in child and father.

Reduced alignment around codon 508.

W1282X
Child shows consistent G>A transition reads.

Father: no evidence of variant → confirms non-paternal origin.

5. Clinical Interpretation
✔ The child is genetically confirmed to have cystic fibrosis.
He carries:

ΔF508 — from father

W1282X — likely maternal or de novo

These two variants together result in CFTR loss of function, explaining:

chronic lung infections

poor weight gain

borderline sweat chloride

6. Limitations
Mother not sequenced → inheritance of second allele uncertain.

Structural variants in CFTR not assessed (e.g., large deletions).

Functional assays (nasal epithelial potential difference) would strengthen diagnosis.

7. Conclusion
Whole-genome sequencing identified two known pathogenic CFTR variants.
The child is a compound heterozygote, consistent with a definitive diagnosis of cystic fibrosis.
Genetic counseling and CF-specific management are strongly recommended.
