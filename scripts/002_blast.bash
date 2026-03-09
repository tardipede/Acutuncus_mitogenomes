#!/usr/bin/bash

# Activate conda nv
conda activate acutuncus_mitogenomes

# Convert the fastq reads files that were produced duting mitogenome assembly into fasta
# The files *_mitomatched.fastq are produced during the run of the script 001_assembly.bash
# They should be moved into a new folder, then move into that folcer with the command cd
seqtk seq -a Acu_gio_RIT094_mitomatched.fastq > Acu_gio_RIT094_mitomatched.fasta
seqtk seq -a Acu_gio_RPL027_mitomatched.fastq > Acu_gio_RPL027_mitomatched.fasta 
seqtk seq -a Acu_gio_RSE010_mitomatched.fastq > Acu_gio_RSE010_mitomatched.fasta 
seqtk seq -a Acu_mec_RSE086_mitomatched.fastq > Acu_mec_RSE086_mitomatched.fasta
seqtk seq -a Acu_mec_RSE087_mitomatched.fastq > Acu_mec_RSE087_mitomatched.fasta

# The *.fasta files used to create the blast databases are the correct mitogenome assemblies obtained by NOVOplasty
# COPY them into the current folder before running the lines below
makeblastdb -in Acu_gio_RIT094.fasta -dbtype nucl -out ./blast_dbs/RIT094
makeblastdb -in Acu_gio_RPL027.fasta -dbtype nucl -out ./blast_dbs/RPL027
makeblastdb -in Acu_gio_RSE010.fasta -dbtype nucl -out ./blast_dbs/RSE010
makeblastdb -in Acu_mec_RSE087.fasta -dbtype nucl -out ./blast_dbs/RSE087
makeblastdb -in Acu_mec_RSE086.fasta -dbtype nucl -out ./blast_dbs/RSE086

# Run the blast searches
blastn -query Acu_gio_RIT094_mitomatched.fasta -db ./blast_dbs/RIT094 -out Acu_gio_RIT094_all_blastn.txt \
       -evalue 1e-5 -num_threads 12 -perc_identity 95 \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

blastn -query Acu_gio_RPL027_mitomatched.fasta -db ./blast_dbs/RPL027 -out Acu_gio_RPL027_all_blastn.txt \
       -evalue 1e-5 -num_threads 12 -perc_identity 95 \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

blastn -query Acu_gio_RSE010_mitomatched.fasta -db ./blast_dbs/RSE010 -out Acu_gio_RSE010_all_blastn.txt \
       -evalue 1e-5 -num_threads 12 -perc_identity 95 \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

blastn -query Acu_mec_RSE087_mitomatched.fasta -db ./blast_dbs/RSE087 -out Acu_mes_RSE087_all_blastn.txt \
       -evalue 1e-5 -num_threads 12 -perc_identity 95 \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

blastn -query Acu_mec_RSE086_mitomatched.fasta -db ./blast_dbs/RSE087 -out Acu_mes_RSE086_all_blastn.txt \
       -evalue 1e-5 -num_threads 12 -perc_identity 95 \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
	   
# Remove the fasta files to free space
rm *.fasta
seqkit stats *.fasta --all > fasta_stats_mitomatched.txt
