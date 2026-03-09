#!/usr/bin/bash

FILENAME=Acu_gio_RIT094.fastq.gz
BASENAME=$(basename $FILENAME |cut -d"." -f1)

# Activate conda env
conda activate acutuncus_mitogenomes

# Decompress reads file (if you start with a gz file, otherwise just skip this line)
gunzip $FILENAME
FILENAME=${BASENAME}.fastq

# Split the input file into 10 parts otherwise porechop will fail
seqkit split2 -p 10 ${FILENAME} -o ${BASENAME}_split
# Run porechopper_abi
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_001.fastq -o ./${FILENAME}.split/${BASENAME}_001_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_002.fastq -o ./${FILENAME}.split/${BASENAME}_002_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_003.fastq -o ./${FILENAME}.split/${BASENAME}_003_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_004.fastq -o ./${FILENAME}.split/${BASENAME}_004_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_005.fastq -o ./${FILENAME}.split/${BASENAME}_005_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_006.fastq -o ./${FILENAME}.split/${BASENAME}_006_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_007.fastq -o ./${FILENAME}.split/${BASENAME}_007_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_008.fastq -o ./${FILENAME}.split/${BASENAME}_008_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_009.fastq -o ./${FILENAME}.split/${BASENAME}_009_chopped.fastq --threads $CPU --check_reads 1000
porechop_abi -i ./${FILENAME}.split/${BASENAME}_split.part_010.fastq -o ./${FILENAME}.split/${BASENAME}_010_chopped.fastq --threads $CPU --check_reads 1000
cat ./${FILENAME}.split/*_chopped.fastq > ./${BASENAME}_chopped.fastq

# Remove non chopped and intermediate files to save space
rm ${FILENAME}
rm -r ./${FILENAME}.split

