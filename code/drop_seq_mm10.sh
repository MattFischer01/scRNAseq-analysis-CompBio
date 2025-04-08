#!/bin/bash

#download STAR to /home/project6/tools (tools folder)
#Example command to run: nohup /home/2025/mfischer10/scRNAseq-analysis-CompBio/code/drop_seq_p2.sh -i /home/2025/mfischer10/scRNAseq-analysis-CompBio/code/reads-to-qc.txt  -o align-quant-output &

if [[ "$1" == "-i" || "$1" == "--input" ]]; then
    input_file="$2"
else
    echo "Unknown parameter passed: $1"
    exit 1
fi

if [[ "$3" == "-o" || "$3" == "--output" ]]; then
    output_dir="$4"
else
    echo "Unknown parameter passed: $3"
    exit 1
fi

if [[ "$5" == "-t" || "$3" == "--tools" ]]; then
    tools_dir="$6"
else
    echo "Unknown parameter passed: $3"
    exit 1
fi

STAR=$(find ${tools_dir} -maxdepth 1 -name "*STAR*" -exec basename {} \;)
DropSeq=$(find ${tools_dir} -maxdepth 1 -name "*dropseq*" -exec basename {} \;)
Picard=$(find ${tools_dir} -maxdepth 1 -name "*picard*" -exec basename {} \;)

STAR_path="${tools_dir}/${STAR}/bin/Linux_x86_64/STAR"
DropSeq_path="${tools_dir}/${DropSeq}"
Picard_path="${tools_dir}/${Picard}"

while IFS=$'\t' read -r read1 read2 output group number; do
    qc_output="${output}/${group}/${group}_${number}"
    align_out_folder="../${output_dir}/${group}/${group}_${number}"
    log_output="../${output_dir}/${group}/${group}_${number}/log"
    log_file="../${output_dir}/${group}_${number}_log.txt"
    read_in="${qc_output}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.fastq"
    read_bam="${qc_output}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.bam"
    out_loc="${align_out_folder}/${group}_${number}"
    group_num="${group}_${number}"

    mkdir -p "${align_out_folder}" "$log_output"

    echo "Processing: $read1 and $read2 > alignment output: $align_output, Logs: $log_output"

#create a genome index for STAR
#we need to start a new directory to hold the indices

#overhang = $(cat max_read.txt)

#if ! test -d "home/project6/genome_index"; then
  #mkdir genome_index

  #/home/project6/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
  #--runThreadN 2 \
  #--runMode genomeGenerate \
  #--genomeDir genome_index \
  #--genomeFastaFiles sample_data/GCF_000001635.26_GRCm38.p6_genomic.fna \
  #--sjdbGTFfile sample_data/GCF_000001635.26_GRCm38.p6_genomic.gtf \
  #--sjdbOverhang 79
#fi

#let's perform an alignment using STAR
$STAR_path \
--genomeDir /home/project6/genome_index \
--readFilesIn $read_in \
--outFileNamePrefix "${out_loc}_star"

#sort output of alignment in queryname order
java -Xmx4g -jar $Picard_path SortSam \
    I="${out_loc}_starAligned.out.sam" \
    O="${out_loc}_aligned.sorted.bam" \
    SO=queryname

java -jar $Picard_path CreateSequenceDictionary \
    R=../sample_data/mm10/mm10.fasta \
    O=../sample_data/mm10/mm10.dict

#merge alignment BAM with unaligned
java -Xmx4g -jar $Picard_path MergeBamAlignment \
    REFERENCE_SEQUENCE=../sample_data/mm10/mm10.fasta \
    UNMAPPED_BAM=$read_bam \
    ALIGNED_BAM="${out_loc}_aligned.sorted.bam" \
    OUTPUT="${out_loc}_merged.bam" \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false

#newer function: tags overlaps with exons and introns - better mapping of exon-exon junctions
"${DropSeq_path}/TagReadWithGeneFunction" \
    I="${out_loc}_merged.bam" \
    O="${out_loc}_star_gene_exon_tagged.bam" \
    ANNOTATIONS_FILE=../sample_data/mm10/mm10.gtf

mkdir "${align_out_folder}/bead_errors"

#detect bead synthesis errors
"${DropSeq_path}/DetectBeadSubstitutionErrors" \
    I="${out_loc}_star_gene_exon_tagged.bam" \
    O="${align_out_folder}/bead_errors/${group_num}_my_clean_subtitution.bam" \
    OUTPUT_REPORT="${align_out_folder}/bead_errors/${group_num}_my_clean.substitution_report.txt"

#detect bead synthesis errors: need to add primer sequences to this
"${DropSeq_path}/DetectBeadSynthesisErrors" \
    I="${align_out_folder}/bead_errors/${group_num}_my_clean_subtitution.bam" \
    O="${align_out_folder}/bead_errors/${group_num}_my_clean.bam" \
    REPORT="${align_out_folder}/bead_errors/${group_num}_my_clean.indel_report.txt" \
    OUTPUT_STATS="${align_out_folder}/bead_errors/${group_num}_my.synthesis_stats.txt" \
    SUMMARY="${align_out_folder}/bead_errors/${group_num}_my.synthesis_stats.summary.txt"

"${DropSeq_path}/BamTagHistogram" \
    I="${out_loc}_star_gene_exon_tagged.bam" \
    O="${out_loc}_out_cell_readcounts.txt.gz" \
    TAG=XC

gunzip "${out_loc}_out_cell_readcounts.txt.gz"

wd=$(pwd)

python3 cell_barcodes.py -i "${out_loc}_out_cell_readcounts.txt" -o  "${wd}/${out_loc}"

#use DigitalExpression to get gene expression matrix
"${DropSeq_path}/DigitalExpression" \
    I="${out_loc}_star_gene_exon_tagged.bam" \
    O="${out_loc}_out_gene_exon_tagged.dge.txt.gz" \
    SUMMARY="${out_loc}_out_gene_exon_tagged.dge.summary.txt" \
    NUM_CORE_BARCODES=100 \
    MIN_NUM_GENES_PER_CELL=500 \
    CELL_BC_FILE="${out_loc}_barcodes.txt"

gunzip "${out_loc}_out_gene_exon_tagged.dge.txt.gz"

done < "$input_file"