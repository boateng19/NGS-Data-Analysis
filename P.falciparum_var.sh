#!/bin/bash
#
# this is a command log to analyse sanger sequence data
# Adapted from:
# https://www.gear-genomics.com/docs/tracy/cli/#basecalling-a-chromatogram-trace-file
#
#
mkdir variant_calls
#
bgzip PF3D7_1460900.fasta
ref=PF3D7_1460900.fasta.gz
#bgzip $ref
#
tracy index $ref
samtools faidx $ref
#
#
#
####_4
######## Variant Calling
######## and convert bcf to vcf file
#
for i in `ls *.ab1 | sed -e 's/\.ab1//'| sort | uniq`;
do
echo ${i}
tracy decompose -v -a plasmodium_falciparum -g $ref -f align -o ./variant_calls/${i}.bcf ${i}.ab1;
done
#
#
########## normalizing bcf files
#
cd variant_calls
#
ref=../PF3D7_1460900.fasta.gz
#
mkdir vcf_files
#
#
for i in `ls *.bcf | sed -e 's/\.bcf.bcf//' | sort | uniq`;
do
echo $i
bcftools norm -O b -o vcf_files/${i}.norm.bcf -f $ref ${i}.bcf.bcf;
done
#
#
#
######### converting normalized bcfs to vcf files
#
cd vcf_files
#
for i in `ls *norm.bcf`;
do
bcftools index ${i};
done
#
for i in `ls *norm.bcf | sed -e 's/\.norm.bcf//' | sort | uniq`;
do
echo $i
bcftools convert -O v -o ${i}.vcf ${i}.norm.bcf
done
#
#
#
######## quality checking vcf files
######## deleting size zero vcf files
######## and removing low qaulity alleles
######## and keeping snps and inDels
#
mkdir fasta_files
#
ref=../../PF3D7_1460900.fasta.gz
#
## delete empty .vcf files
find . -size 0 -print -delete
#
for i in `ls *.vcf | sed -e 's/\.vcf//' | sort | uniq`;
do
bgzip ${i}.vcf
tabix ${i}.vcf.gz
# remove low quality alleles
bcftools view -i '%QUAL>=20' ${i}.vcf.gz > ${i}_filter.vcf
bgzip ${i}_filter.vcf
tabix ${i}_filter.vcf.gz

# keep snvs and inDels
bcftools view --types snps ${i}_filter.vcf.gz -Ov -o ${i}_filter_final.vcf

bgzip ${i}_filter_final.vcf
tabix ${i}_filter_final.vcf.gz

#combining vcf files into one
cat ${i}_filter_final.vcf.gz > all_data.vcf

############# creating consensus fasta
#
samtools faidx $ref
#
bcftools index ${i}_filter_final.vcf.gz
bcftools consensus -f $ref ${i}_filter_final.vcf.gz > ./fasta_files/${i}.fasta
done
cat ./fasta_files/*.fasta > ./fasta_files/all_arps.fasta
#
#
#
#
########## changing chromosome name
#
echo "PF3D7_1460900 NC_037283.1" > newname.txt
bcftools annotate --rename-chrs newname.txt all_data.vcf -o all_data_renamed.vcf
#
#
############# snpEff annotation
mkdir varANN
#
#
java -Xmx8g -jar /home/prince/snpEff/snpEff.jar ann -v NC_037283 all_data_renamed.vcf > ./varANN/all_data_ann.vcf
#
cd varANN
#
cat all_data_ann.vcf | perl /home/prince/snpEff/scripts/vcfEffOnePerLine.pl |\
java -Xmx8g -jar /home/prince/snpEff/SnpSift.jar extractFields -e "NA" - \
CHROM POS REF ALT QUAL FILTER "ANN[*].BIOTYPE" "ANN[*].EFFECT" "ANN[*].HGVS_P" "ANN[*].ERRORS" > all_data_ANN.txt
#



