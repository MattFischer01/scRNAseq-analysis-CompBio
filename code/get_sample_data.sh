#Code to extract the first 10000 lines of naive SRR data so I can troubleshoot QC before proceeding with all data. 
#Naive: SRR8878106_1
zcat /home/project6/scRNA_SRA/naive/SRR8878106_1.fastq.gz | head -n 40000 | gzip > /home/project6/sample_data/SRR8878106_1.fastq.gz
zcat /home/project6/scRNA_SRA/naive/SRR8878106_2.fastq.gz | head -n 40000 | gzip > /home/project6/sample_data/SRR8878106_2.fastq.gz

#Dropseq sequencing libraries produce paired-end reads: read 1 contains both a cell barcode and a molecular barcode (UMI); read 2 is aligned to the reference genome.

#check if the paired fastq files are in the same format
zcat /home/project6/sample_data/SRR8878106_1.fastq.gz | head -n 8
# @SRR8878106.1 1 length=20
# AATTCGTCCCTATCGCAAAA
# +SRR8878106.1 1 length=20
# FFFFF,FF,F,FFF:,FF:F
# @SRR8878106.2 2 length=20
# CTATGCCCTAACCGCTCATG
# +SRR8878106.2 2 length=20
# FFFFFFFFFFFFFFF,FFFF

zcat /home/project6/sample_data/SRR8878106_2.fastq.gz | head -n 8
# @SRR8878106.1 1 length=80
# GCTATAAGAAGGGGATAGAAACAAACAACAGATGGCTAGGGTCAGAACTCAATGGTTGAGCACTTACCTCCTATATGTAA
# +SRR8878106.1 1 length=80
# FFFFFFFFF,FFFF:F:,FFFFFFFFFFF,,FFFFFF:,FFFFFF:FF::FFF,:FFFFFFFFF:FF:F,F,F,FF:FFF
# @SRR8878106.2 2 length=80
# CTGCTGCATTGCTGCCCTGCACCAAACATGCCTAGGCCGACGAGTTCCCAGTTAAGTCGTATAACCTGGCTCCAGTGTGT
# +SRR8878106.2 2 length=80
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF

#Looks about right I believe
#he first 12 bases of the barcoded read contain the cell barcode, so we’ll copy those bases over to a BAM tag (XC) on the genome read. Then we’ll take the next 8 bases containing the molecular barcode and copy them over as another BAM tag (XM). 

#Rename files for ease:
mv SRR8878106_1.fastq.gz naive_sample_01.fastq.gz
mv SRR8878106_2.fastq.gz naive_sample_02.fastq.gz



#Now I need to convert the fastq files into a Picard-queryname-sorted BAM file. I will use Picard FastqToSam:
#https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam

#The code continues with drop_seq_qc.sh
