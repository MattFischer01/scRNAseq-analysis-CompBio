#!/bin/bash

#This script will allow the user to download dependencies necessary to run the pipeline. They should run this script before the pipeline, so they will only have to do it one. The tools will be downloaded in a tool directory found in the parent directory. 

#This code should be run in ./code directory. 

#Make directory for tools. 
mkdir -p ../tools

#First tool: Picard
#Download picard in tools directory 
wget -P ../tools https://github.com/broadinstitute/picard/releases/latest/download/picard.jar

#Second tool: Drop-seq
#Download drop-seq in tools directory 
wget -P ../tools https://github.com/broadinstitute/Drop-seq/releases/download/v3.0.2/dropseq-3.0.2.zip

unzip ../tools/dropseq-3.0.2.zip -d ../tools

#Third tool: STAR
#download STAR in the tools directory
wget -P ../tools https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
cd ../tools
tar -xzf ./2.7.11b.tar.gz

#compile STAR
cd ./STAR-2.7.11b/source
make STAR

#cd back into the code directory
cd ../../../code

#Download necessary metadata to run STAR alignment.
wget -P ../meta_data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz 
wget -P ../meta_data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz 
mv ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fna.gz ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta.gz

gunzip ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta.gz
gunzip ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz
