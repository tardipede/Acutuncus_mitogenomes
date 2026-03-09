#!/usr/bin/bash

FILENAME=Acu_gio_RIT094.fastq.gz
BASENAME=$(basename $FILENAME |cut -d"." -f1)

# Activate conda env
conda activate acutuncus_mitogenomes

# Decompress reads file (if you start with a gz file, otherwise just skip this line)
gunzip $FILENAME
FILENAME=${BASENAME}.fastq

# create a diamond-formatted database file
diamond makedb --in ../databases/tardigrada_mito_proteins.db.v01.fas -d ../databases/diamond_db

# running a search in blastx mode
diamond blastx -d ../databases/diamond_db --query-gencode 5 --sensitive --max-target-seqs 1 --threads 12 -q ${FILENAME} --outfmt 6 qseqid -o matches.tsv

# Extract matched sequences and filter by length
seqtk subseq ${FILENAME} matches.tsv > ${BASENAME}_mitomatched.fastq

##### SEPARATE INTO ILLUMINA-LIKE fragments
seqkit sliding ${BASENAME}_mitomatched.fastq -s 150 -W 150 -o ${BASENAME}_sliding.fastq
seqkit seq ${BASENAME}_sliding.fastq -Q 30 > illuminalike.fastq

# Assembly with novoplasty
NOVOPlasty4.3.5.pl -c config.txt

# Assembly with flye standard mode
flye --nano-raw ${BASENAME}_mitomatched.fastq -t 10 -i 1 -o ./flye_assembly

# Assembly with flye meta mode
flye --nano-raw ${BASENAME}_mitomatched.fastq --meta -t 10 -i 1 -o ./flye_assembly_meta

# Assembly with flye with genome size specified
flye --nano-raw ${BASENAME}_mitomatched.fastq -g 14500 --asm-coverage 1000 -t 10 -i 1 -o ./flye_assembly_gsize

# Palindrome detection with alignment size 1000
# Code from https://doi.org/10.1093/nar/gkad647
CPU=12

# Create new folder where to store the palindrome treated reads
mkdir palindrome
cd palindrome

# Mode the mitomatched reads to the palindrome folder
mv ../${BASENAME}_mitomatched.fastq ./${BASENAME}_mitomatched.fastq

# Run the palindrome finding script with alignment length of 1000
minimap2 -t $CPU -x ava-ont ./${BASENAME}_mitomatched.fastq ./${BASENAME}_mitomatched.fastq | cat |../../softwares/palindrome_detection/pafIdentifyPalimdrom.pl 1> palimProp.1stIte.list 2> 1stIte.log
../../softwares/palindrome_detection/fastq_partition.and.chop.palindrome.pl palimProp.1stIte.list ./${BASENAME}_mitomatched.fastq 1000
minimap2 -t $CPU -x ava-ont ${BASENAME}_mitomatched.fastq.exclude.fq.gz ${BASENAME}_mitomatched.fastq.exclude.fq.gz | cat | ../../softwares/palindrome_detection/pafIdentifyPalimdrom.pl 1> palimProp.2ndIte.list 2> 2ndIte.log
../../softwares/palindrome_detection/fastq_partition.and.chop.palindrome.pl palimProp.2ndIte.list ${BASENAME}_mitomatched.exclude.fq.gz 1000
cat ${BASENAME}_mitomatched.fastq.include.fq.gz  ${BASENAME}_mitomatched.exclude.fq.gz.include.fq.gz > ${BASENAME}_mitomatched.palindrome_treated1000.fq.gz


# Run the palindrome finding script with alignment length of 100
# The firts itaration of minimap2 is skipped as results from the previous rrun are used to save time
../../softwares/palindrome_detection/fastq_partition.and.chop.palindrome.pl palimProp.1stIte.list ./${BASENAME}_mitomatched.fastq 100
minimap2 -t $CPU -x ava-ont ${BASENAME}_mitomatched.fastq.exclude.fq.gz ${BASENAME}_mitomatched.fastq.exclude.fq.gz | cat | ../../softwares/palindrome_detection/pafIdentifyPalimdrom.pl 1> palimProp.2ndIte.list 2> 2ndIte.log
../../softwares/palindrome_detection/fastq_partition.and.chop.palindrome.pl palimProp.2ndIte.list ${BASENAME}_mitomatched.exclude.fq.gz 100
cat ${BASENAME}_mitomatched.fastq.include.fq.gz  ${BASENAME}_mitomatched.exclude.fq.gz.include.fq.gz > ${BASENAME}_mitomatched.palindrome_treated100.fq.gz

# Remove heavy intermediate files
rm ${BASENAME}_mitomatched.fastq.exclude.fq.gz
rm ${BASENAME}_mitomatched.fastq.include.fq.gz
rm ${BASENAME}_mitomatched.exclude.fq.gz.include.fq.gz
rm ${BASENAME}_mitomatched.exclude.fq.gz.exclude.fq.gz
rm palimProp.2ndIte.list
rm palimProp.1stIte.list
rm 2ndIte.log
rm 1stIte.log



# Assembly with flye standard mode
flye --nano-raw ${BASENAME}_mitomatched.palindrome_treated1000.fq.gz -t 10 -i 1 -o ./flye_assembly_pal1000

# Assembly with flye meta mode
flye --nano-raw ${BASENAME}_mitomatched.palindrome_treated1000.fq.gz --meta -t 10 -i 1 -o ./flye_assembly_meta_pal1000

# Assembly with flye with genome size specified
flye --nano-raw ${BASENAME}_mitomatched.palindrome_treated1000.fq.gz -g 14500 --asm-coverage 1000 -t 10 -i 1 -o ./flye_assembly_gsize_pal1000

# Assembly with flye standard mode
flye --nano-raw ${BASENAME}_mitomatched.palindrome_treated100.fq.gz -t 10 -i 1 -o ./flye_assembly_pal100

# Assembly with flye meta mode
flye --nano-raw ${BASENAME}_mitomatched.palindrome_treated100.fq.gz --meta -t 10 -i 1 -o ./flye_assembly_meta_pal100

# Assembly with flye with genome size specified
flye --nano-raw ${BASENAME}_mitomatched.palindrome_treated100.fq.gz -g 14500 --asm-coverage 1000 -t 10 -i 1 -o ./flye_assembly_gsize_pal100

