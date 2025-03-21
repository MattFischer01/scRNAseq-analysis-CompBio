#Step 1: Now that I have the necessary metadata. I need to convert the fastq files into a Picard-queryname-sorted BAM file. I will use Picard FastqToSam:
#https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam

#download Picard: wget https://github.com/broadinstitute/picard/releases/latest/download/picard.jar

java -jar /home/project6/tools/picard.jar FastqToSam \
       F1=/home/project6/sample_data/naive_sample_01.fastq.gz \
       F2=/home/project6/sample_data/naive_sample_02.fastq.gz \
       O=/home/project6/dropseq-qc-output/naive_sample_unaligned_read_pairs.bam \
       SM=naive001 \
       RG=naive

#Step 2: TagBamWithReadSequenceExtended. This code extracts the cell and molecular barcode encoding read and creates a new BAM tag with those bases on the genome read. XM tag for molecular barcodes. XC for cell barcodes. I need to run the program once per barcode extraction to add a tag, first the cell barcode to extract bases 1-12, then run again to extract the molecular barcode fomr bases 13-20 of the barcode read. This second run has an option to drop the barcode read (DISCARD_READ), which makes the output BAM have unpaired reads with additional tags. 

#Step 2a: Cell barcode
/home/project6/tools/dropseq-3.0.2/TagBamWithReadSequenceExtended \
       INPUT=/home/project6/dropseq-qc-output/naive_sample_unaligned_read_pairs.bam \
       OUTPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_Cell.bam \
       SUMMARY=/home/project6/dropseq-qc-output/naive_unaligned_tagged_Cellular.bam_summary.txt \
       BASE_RANGE=1-12 \
       BASE_QUALITY=10 \
       BARCODED_READ=1 \
       TAG_NAME=XC \
       NUM_BASES_BELOW_QUALITY=1

#Step 2b: Molcular barcode
/home/project6/tools/dropseq-3.0.2/TagBamWithReadSequenceExtended \
       INPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_Cell.bam \
       OUTPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_CellMolecular.bam \
       SUMMARY=/home/project6/dropseq-qc-output/naive_unaligned_tagged_Molecular.bam_summary.txt \
       BASE_RANGE=13-20 \
       BASE_QUALITY=10 \
       BARCODED_READ=1 \
       DISCARD_READ=TRUE \
       TAG_NAME=XM \
       NUM_BASES_BELOW_QUALITY=1

#Step 3: FilterBam:
##Filtering step to remove read where the cell or molecular barcode has low quality bases (based on the XQ assigned tag during TagBamWithReadSequenceExtended)
/home/project6/tools/dropseq-3.0.2/FilterBam \
       TAG_REJECT=XQ \
       INPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_CellMolecular.bam \
       OUTPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_filtered.bam

#Results: Total 10000 reads processed.  10000 reads accepted; 0 reads rejected.

#Step 4: TrimStartingSequence:
##One of two sequence cleanup programs designed to rim away any extra sequence that might have snuck it's way into the reads. In this case, we trim the SMART adapter that can occur 5' of the read. 

/home/project6/tools/dropseq-3.0.2/TrimStartingSequence \
       INPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_filtered.bam \
       OUTPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_trimmed_smart.bam \
       OUTPUT_SUMMARY=/home/project6/dropseq-qc-output/naive_adapter_trimming_report.txt \
       SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
       MISMATCHES=0 \
       NUM_BASES=5

#Results:  Number of reads trimmed: 132 total reads: 10000

#Step 5: PolyATrimmer 
##Second sequence cleanup program designed to trim away trailing polyA tails from reads:
/home/project6/tools/dropseq-3.0.2/PolyATrimmer \
       INPUT=/home/project6/dropseq-qc-output/naive_unaligned_tagged_trimmed_smart.bam \
       OUTPUT=/home/project6/dropseq-qc-output/naive_unaligned_mc_tagged_polyA_filtered.bam \
       OUTPUT_SUMMARY=/home/project6/dropseq-qc-output/naive_polyA_trimming_report.txt \
       MISMATCHES=0 \
       NUM_BASES=6 

#Results
# INFO    2025-03-20 10:55:18     PolyATrimmer    Total 10000 reads processed.
# INFO    2025-03-20 10:55:18     PolyATrimmer    Number of reads trimmed: 960
# INFO    2025-03-20 10:55:18     PolyATrimmer    Number of reads completely trimmed: 11

#Step 6: SamToFastq
##Now that we have qced the data, it is time to align, we first extract the FASTQ files using Picards SamToFastq program:
java -Xmx4g -jar /home/project6/tools/picard.jar SamToFastq \
       INPUT=/home/project6/dropseq-qc-output/naive_unaligned_mc_tagged_polyA_filtered.bam \
       FASTQ=/home/project6/dropseq-qc-output/naive_unaligned_mc_tagged_polyA_filtered.fastq


#Step 7: Alignment to STAR - Elise
