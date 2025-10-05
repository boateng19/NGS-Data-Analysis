#!/bin/bash

# Differential Expression Analysis for UV Stress in Arabidopsis Vasculature
# Berkowitz et al. (2021) dataset

# Create project directory structure
mkdir -p vasculature_uv_analysis/{raw_data,processed_data,results}
cd vasculature_uv_analysis

# Define sample IDs
control_samples=("SRR12808527" "SRR12808528" "SRR12808529")
uv_samples=("SRR12808497" "SRR12808498" "SRR12808499")

# Download SRA files using fasterq-dump
echo "Downloading control samples..."
for sample in "${control_samples[@]}"; do
    fasterq-dump --split-files $sample -O raw_data/
done

echo "Downloading UV-treated samples..."
for sample in "${uv_samples[@]}"; do
    fasterq-dump --split-files $sample -O raw_data/
done

# Quality control with FastQC
echo "Running FastQC..."
mkdir -p results/fastqc
fastqc raw_data/*.fastq -o results/fastqc/

# Trimming with Trimmomatic
echo "Trimming adapters and low-quality bases..."

mkdir -p processed_data/trimmed reports

for fq in raw_data/*.fastq.gz; do
  base=$(basename "$fq" .fastq.gz)
  fastp \
    -i "$fq" \
    -o processed_data/trimmed/"${base}.trim.fastq.gz" \
    --html reports/"${base}_report.html" \
    --json reports/"${base}_report.json"
done

# Download Arabidopsis thaliana reference genome
echo "Downloading reference genome..."
mkdir -p reference
wget -P reference/ https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget -P reference/ https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-52/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gtf.gz

gunzip reference/*.gz

# Build HISAT2 index
echo "Building HISAT2 index..."
#hisat2-build reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa reference/athaliana_index
hisat2-build Arabidopsis_thaliana.TAIR10.dna.toplevel.fa tair10_index

# Alignment with HISAT2
echo "Aligning reads to reference genome..."
mkdir -p processed_data/aligned

for fq in processed_data/trimmed/*.trim.fastq.gz; do
  base=$(basename "$fq" .trim.fastq.gz)
  echo "Aligning $base ..."
  hisat2 -x reference/tair10_index \
    -U "$fq" \
    -S processed_data/aligned/${base}.sam \
    --threads 4

done
# Convert SAM to BAM and sort

echo "Converting and sorting BAM files"
for fq in processed_data/aligned/*.sam; do
  base=$(basename "$fq" .sam)   
  echo "Processing $base ..."

  samtools view -@ 4 -bS "$fq" \
    | samtools sort -@ 4 -o processed_data/aligned/${base}.sorted.bam

  samtools index processed_data/aligned/${base}.sorted.bam

done


# Count reads with featureCounts
echo "Counting reads per gene..."
mkdir -p counts

featureCounts -T 4 -s 0 \
  -a reference/Arabidopsis_thaliana.TAIR10.52.gtf \
  -o counts/gene_counts.txt \
  processed_data/aligned/*.sorted.bam

echo "Preprocessing complete! Proceed to R analysis."


# Prepare count matrix for DESeq2

cut -f1,7- results/gene_counts.txt > counts/counts_matrix.txt

