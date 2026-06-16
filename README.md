# Acutuncus_mitogenomes
Assembly of Acutuncus mitogenomes with a quick solution to deal with MDA generated palindromes  
Associated with the article:  
Vecchi, Surmacz, Jonsson & Stec (2026) **Fragmentation of Long Reads Enables Reliable Mitogenome Assembly From Whole-Genome Amplification Data With Pervasive Palindromic Reads.** *Molecular Ecology Resources*, 26:e70165. https://doi.org/10.1111/1755-0998.70165

## 1) Create a folder name data and place there the raw reads files
Download the raw reads file from the NCBI project PRJNA1287536.
The runs names should be downloaded with their original names, otherwise rename them as below:  
    Acu_gio_RIT094.fastq  
    Acu_gio_RPL027.fastq  
    Acu_gio_RSE010.fastq  
    Acu_mec_RSE086.fastq  
    Acu_mec_RSE087.fastq  

## 2) Create conda env
```
conda env create -f ./acutuncus_mitogenomes.yml
```

## 3) Run scripts in order
    * scripts/001_assembly.bash (Run in bash - must be separately run for each individual).
    * scripts/002_blast.bash (Run in bash - before running it check how to prepare the data as written inside the script file).
    * scripts/003_mitos.bash (Run in bash - before running it check how to prepare the data as written inside the script file).
    * scripts/004_blast_results_analysis.r (Run in R - before running it check how to prepare the data as written inside the script file).

## 4) Notes before running script 001_assembly
The config and COI reference files needed to run NOVOplasty are provided in the folder *NOVOplasty_config_material*, before running the script 001 place in the working directory the appropriate reference sequence in fasta format and the config file. Inside the config file change the *Seed Input* field to match the name of the reference fasta file used. If you only want to assembly the mitogenome with NOVOplasty, stop at line 27 of the script, from line 28 onward it does the 9 different assemblies with Flye.

## 5) Reads primers and adapters chopping
The reads present in NCBI have been already cleaned from adapters and barcodes, however the code used is available to check at *scripts/999_chopper.bash*.
