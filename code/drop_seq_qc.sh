#!/bin/bash

#This code will read in a text file and run paired end reads through Drop-seq's qc process.

#Example command to run: nohup /home/2025/mfischer10/scRNAseq-analysis-CompBio/code/drop_seq_qc.sh -i /home/2025/mfischer10/scRNAseq-analysis-CompBio/code/reads-to-qc.txt &

#Takes an argument for the path to the text file and assigns it to variable input_file.
if [[ "$1" == "-i" || "$1" == "--input" ]]; then
    input_file="$2"
else
    echo "Unknown parameter passed: $1"
    exit 1
fi


while IFS=$'\t' read -r read1 read2 output group number; do
       qc_output="${output}/${group}/${group}_${number}"
       log_output="${output}/${group}/${group}_${number}/log"
       log_file="${log_output}/${group}_${number}_log.txt"

       mkdir -p "$qc_output" "$log_output"

       echo "Processing: $read1 and $read2 > QC output: $qc_output, Logs: $log_output"


#Step 1: I need to convert the fastq files into a Picard-queryname-sorted BAM file. I will use Picard FastqToSam:
#https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam

#download Picard: wget https://github.com/broadinstitute/picard/releases/latest/download/picard.jar

{
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Picard FastqToSam"
       java -jar /home/project6/tools/picard.jar FastqToSam \
              F1="$read1" \
              F2="$read2" \
              O="$qc_output/${group}_${number}_unaligned_read_pairs.bam" \
              SM="${group}_${number}" \
              RG="${group}"
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Picard FastqToSam"
} >> "$log_file" 2>&1


#Step 2: TagBamWithReadSequenceExtended. This code extracts the cell and molecular barcode encoding read and creates a new BAM tag with those bases on the genome read. XM tag for molecular barcodes. XC for cell barcodes. I need to run the program once per barcode extraction to add a tag, first the cell barcode to extract bases 1-12, then run again to extract the molecular barcode fomr bases 13-20 of the barcode read. This second run has an option to drop the barcode read (DISCARD_READ), which makes the output BAM have unpaired reads with additional tags. 

#Step 2a: Cell barcode
{
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Dropseq TagBamWithReadSequenceExtended - Cell Barcode"
       /home/project6/tools/dropseq-3.0.2/TagBamWithReadSequenceExtended \
              INPUT="$qc_output/${group}_${number}_unaligned_read_pairs.bam" \
              OUTPUT="$qc_output/${group}_${number}_unaligned_tagged_Cell.bam" \
              SUMMARY="$qc_output/${group}_${number}_unaligned_tagged_Cellular.bam_summary.txt" \
              BASE_RANGE=1-12 \
              BASE_QUALITY=10 \
              BARCODED_READ=1 \
              TAG_NAME=XC \
              NUM_BASES_BELOW_QUALITY=1
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Dropseq TagBamWithReadSequenceExtended - Cell Barcode"
} >> "$log_file" 2>&1

#Step 2b: Molcular barcode
{
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Dropseq TagBamWithReadSequenceExtended - Molecular Barcode"
       /home/project6/tools/dropseq-3.0.2/TagBamWithReadSequenceExtended \
              INPUT="$qc_output/${group}_${number}_unaligned_tagged_Cell.bam" \
              OUTPUT="$qc_output/${group}_${number}_unaligned_tagged_CellMolecular.bam" \
              SUMMARY="$qc_output/${group}_${number}_unaligned_tagged_Molecular.bam_summary.txt" \
              BASE_RANGE=13-20 \
              BASE_QUALITY=10 \
              BARCODED_READ=1 \
              DISCARD_READ=TRUE \
              TAG_NAME=XM \
              NUM_BASES_BELOW_QUALITY=1
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Dropseq TagBamWithReadSequenceExtended - Molecular Barcode"
} >> "$log_file" 2>&1


#Step 3: FilterBam:
##Filtering step to remove read where the cell or molecular barcode has low quality bases (based on the XQ assigned tag during TagBamWithReadSequenceExtended)
{
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Dropseq FilterBam"
       /home/project6/tools/dropseq-3.0.2/FilterBam \
              TAG_REJECT=XQ \
              INPUT="$qc_output/${group}_${number}_unaligned_tagged_CellMolecular.bam" \
              OUTPUT="$qc_output/${group}_${number}_unaligned_tagged_filtered.bam"
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Dropseq FilterBam"
} >> "$log_file" 2>&1



#Step 4: TrimStartingSequence:
##One of two sequence cleanup programs designed to rim away any extra sequence that might have snuck it's way into the reads. In this case, we trim the SMART adapter that can occur 5' of the read. 

{
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Dropseq TrimStartingSequence"
       /home/project6/tools/dropseq-3.0.2/TrimStartingSequence \
              INPUT="$qc_output/${group}_${number}_unaligned_tagged_filtered.bam" \
              OUTPUT="$qc_output/${group}_${number}_unaligned_tagged_trimmed_smart.bam" \
              OUTPUT_SUMMARY="$qc_output/${group}_${number}_adapter_trimming_report.txt" \
              SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
              MISMATCHES=0 \
              NUM_BASES=5
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Dropseq TrimStartingSequence"
} >> "$log_file" 2>&1


#Step 5: PolyATrimmer 
##Second sequence cleanup program designed to trim away trailing polyA tails from reads:
{
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Dropseq PolyATrimmer"
       /home/project6/tools/dropseq-3.0.2/PolyATrimmer \
              INPUT="$qc_output/${group}_${number}_unaligned_tagged_trimmed_smart.bam" \
              OUTPUT="$qc_output/${group}_${number}_unaligned_mc_tagged_polyA_filtered.bam" \
              OUTPUT_SUMMARY="$qc_output/${group}_${number}_polyA_trimming_report.txt" \
              MISMATCHES=0 \
              NUM_BASES=6 
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Dropseq PolyATrimmer"
} >> "$log_file" 2>&1


#Step 6: SamToFastq
##Now that we have qced the data, it is time to align, we first extract the FASTQ files using Picards SamToFastq program:

{
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Picard SamToFastq"
       java -Xmx4g -jar /home/project6/tools/picard.jar SamToFastq \
              INPUT="$qc_output/${group}_${number}_unaligned_mc_tagged_polyA_filtered.bam" \
              FASTQ="$qc_output/${group}_${number}_unaligned_mc_tagged_polyA_filtered.fastq"
       echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Picard SamToFastq"
} >> "$log_file" 2>&1


done < "$input_file"


#Step 7: Continue with Star alignment in drop_Seq_p2.sh
