#!/bin/bash

# Create a new folder for the MITOS2 run
# Inside this folder copy the assemblies from NOVOplasty
# Download the MITOS2 reference database (refseq89m) from here: https://doi.org/10.5281/zenodo.3685310
# One copy of the MITOS2 refseq89m database is alredy provided in the "databases" folder
# And place the unzipped database folder inide the folder for the MITOS2 runs

# Activate conda environment
conda activate acutuncus_mitogenomes

# Run MITOS2 for each of the 5 assemblies
FASTA_FILE=Acu_gio_RIT094.fasta
BASENAME=$(basename $FASTA_FILE |cut -d"." -f1)
mkdir ./mitos_out_${BASENAME}
runmitos.py --code 5 --refdir ./refseq89m/ -r ./ \
            --input $FASTA_FILE --zip --noplots \
            --outdir ./mitos_out_${BASENAME} --finovl 100 --fragovl 75 --locandgloc

FASTA_FILE=Acu_gio_RPL027.fasta
BASENAME=$(basename $FASTA_FILE |cut -d"." -f1)
mkdir ./mitos_out_${BASENAME}
runmitos.py --code 5 --refdir ./refseq89m/ -r ./ \
            --input $FASTA_FILE --zip --noplots \
            --outdir ./mitos_out_${BASENAME} --finovl 100 --fragovl 75 --locandgloc

FASTA_FILE=Acu_gio_RSE010.fasta
BASENAME=$(basename $FASTA_FILE |cut -d"." -f1)
mkdir ./mitos_out_${BASENAME}
runmitos.py --code 5 --refdir ./refseq89m/ -r ./ \
            --input $FASTA_FILE --zip --noplots \
            --outdir ./mitos_out_${BASENAME} --finovl 100 --fragovl 30 --locandgloc

FASTA_FILE=Acu_mec_RSE087.fasta
BASENAME=$(basename $FASTA_FILE |cut -d"." -f1)
mkdir ./mitos_out_${BASENAME}
runmitos.py --code 5 --refdir ./refseq89m/ -r ./ \
            --input $FASTA_FILE --zip --noplots \
            --outdir ./mitos_out_${BASENAME} --finovl 100 --fragovl 30 --locandgloc

FASTA_FILE=Acu_mec_RSE086.fasta
BASENAME=$(basename $FASTA_FILE |cut -d"." -f1)
mkdir ./mitos_out_${BASENAME}
runmitos.py --code 5 --refdir ./refseq89m/ -r ./ \
            --input $FASTA_FILE --zip --noplots \
            --outdir ./mitos_out_${BASENAME} --finovl 100 --fragovl 30 --locandgloc