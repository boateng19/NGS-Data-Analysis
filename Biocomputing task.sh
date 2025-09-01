# Project 1: Bash Basic
#!#/bin/bash

# 1. Printing my name
echo "Agyekum Boateng"

# 2. Create a folder titled AgyekumBoateng
mkdir  AgyekumBoateng

# 3. Create another directory titled biocomputing and change into it
mkdir biocomputing && cd biocomputing

# 4. Downloading the 3 files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
pe.gbk

# 5. Move the .fna file to the folder titled AgyekumBoateng
mv wildtype.fna ../AgyekumBoateng/

# 6. Delete duplicate gbk file
rm -f wildtype.gbk.1

# 7. Confirming if the .fna file is mutant or wild type
grep "tatatata" ../AgyekumBoateng/wildtype.fna
grep "tata" ../AgyekumBoateng/wildtype.fna

# 8. If mutant, print all matching lines into a new file
grep "tatatata" ../AgyekumBoateng/wildtype.fna > mutant_sequences.txt
cat mutant_sequence.txt

# 9. Counting the number of lines (excluding header) in the .gbk file
sed -n '/^ORIGIN/,/^\/\//p' wildtype.gbk | grep -v "ORIGIN" | grep -v "//" | wc -l

# 10. Printing the sequence length of the .gbk file

grep "LOCUS" wildtype.gbk

# 11. Printing the source organism of the .gbk file

grep "SOURCE" wildtype.gbk

# 12. Listing all gene names in the .gbk file

grep "/gene=" wildtype.gbk


# 13. Clear terminal and print all commands used today
clear
history

# 14. Listing files in the two folders
#files in biocomputing
ls

# files in AgyekumBoateng
cd ../AgyekumBoateng && ls






#Project 2


# 1. Activating Conda base environment
conda activate

#2. Create a conda environment named funtools
conda create -n funtools

#3. Activate the funtools environment
conda activate funtools

#4. Install Figlet using conda
conda install figlet

#5. Run figlet <my name>
figlet Agyekum

#6. Install bwa through the bioconda channel

conda install -c bioconda bwa

#7. Run bwa
bwa

#8. Install blast through the bioconda channel
conda install -c bioconda blast

#9. Run blast
blast

#10. Install samtools through the bioconda channel
conda install -c bioconda samtools

#11. Run samtools
samtools

#12. Install bedtools through the bioconda channel
conda install -c bioconda bedtools

#13. Run bedtools
bedtools

#14. Install spades.py through the bioconda channel
conda install -c bioconda::spades

#15. Run spades
spades

#16. Install bcftools through the bioconda channel
conda install -c bioconda bcftools

#17. Run bcftools
bcftools

#18. Install fastp through the bioconda channel

conda install -c bioconda fastp

#19. Run fastp
fastp

#20. Install multiqc through the bioconda channel
conda install -c bioconda multiqc

#21. Run multiqc
multiqc
