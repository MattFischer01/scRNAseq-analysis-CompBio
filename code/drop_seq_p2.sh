#!/bin/bash

if [[ "$1" == "-i" || "$1" == "--input" ]]; then #take in sample metadata file as input_file
    input_file="$2"
else
    echo "Unknown parameter passed: $1"
    exit 1
fi

if [[ "$3" == "-o" || "$3" == "--output" ]]; then #take in output directory name as output_dir
    output_dir="$4"
else
    echo "Unknown parameter passed: $3"
    exit 1
fi

if [[ "$5" == "-t" || "$5" == "--tools" ]]; then #take in tools directory as tools_dir
    tools_dir="$6"
else
    echo "Unknown parameter passed: $5"
    exit 1
fi

if [[ "$7" == "-x" || "$7" == "--index" ]]; then #take in tools directory as tools_dir
    index_dir="$8"
else
    echo "Unknown parameter passed: $7"
    exit 1
fi

if [[ "$9" == "-g" || "$9" == "--min_genes" ]]; then #take in tools directory as tools_dir
    min_genes="${10}"
else
    echo "Unknown parameter passed: $9"
    exit 1
fi

STAR=$(find ${tools_dir} -maxdepth 1 -name "*STAR*" -exec basename {} \;) #search in tools_dir for name containing STAR
DropSeq=$(find ${tools_dir} -maxdepth 1 -name "*dropseq*" -exec basename {} \;) #search in tools_dir for name containing dropseq
Picard=$(find ${tools_dir} -maxdepth 1 -name "*picard*" -exec basename {} \;) #search in tools_dir for name containing picard

STAR_path="${tools_dir}/${STAR}/bin/Linux_x86_64/STAR" #construct full STAR path
DropSeq_path="${tools_dir}/${DropSeq}" #construct full DropSeq path
Picard_path="${tools_dir}/${Picard}" #construct full Picard path

#this code generates a genome index if it does not already exist
if [[ "$index_dir" == "auto" ]]; then
  mkdir ../genome_index

  $STAR_path \
  --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeDir ../genome_index \
  --genomeFastaFiles ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta \
  --sjdbGTFfile ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.gtf \
  --sjdbOverhang 79

  index_loc="../genome_index"
else
  index_loc="$index_dir"
fi

while IFS=$'\t' read -r read1 read2 output group number; do
    qc_output="${output}/${group}/${group}_${number}" #location to pull QC results from
    align_out_folder="../${output_dir}/${group}/${group}_${number}" #folder to put Align-Quant results
    log_output="../${output_dir}/${group}/${group}_${number}/log" #log directory
    log_file="../${output_dir}/${group}_${number}_log.txt" #log file name
    read_in="${qc_output}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.fastq" #file from QC to read into STAR
    read_bam="${qc_output}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.bam" #file from QC to read into SortSam
    out_loc="${align_out_folder}/${group}_${number}" #output location including prefix
    group_num="${group}_${number}" #just prefix
    path_to_DGE="${group} ${number} ${out_loc}_out_gene_exon_tagged.dge.txt"

    mkdir -p "${align_out_folder}" "$log_output" #make the align-quant output directory and log directory

    echo "Processing: $read1 and $read2 > alignment output: $align_out_folder, Logs: $log_output"

#let's perform an alignment using STAR: 
    #note: index is hard-coded for now, but the user should either supply the location to this as input or have the code generate it
$STAR_path \
--genomeDir $index_loc \
--readFilesIn $read_in \
--outFileNamePrefix "${out_loc}_star"

#sort output of alignment in queryname order
java -Xmx4g -jar $Picard_path SortSam \
    I="${out_loc}_starAligned.out.sam" \
    O="${out_loc}_aligned.sorted.bam" \
    SO=queryname \
    TMP_DIR=/home/2025/estagaman/tmp \
    QUIET=true

#####creating metadata from the mouse reference files which can be used later
java -jar $Picard_path CreateSequenceDictionary \
    R=../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta \
    O=../meta_data/GCF_000001635.26_GRCm38.p6_genomic.dict \
    QUIET=true

#merge alignment BAM with unaligned
java -Xmx4g -jar $Picard_path MergeBamAlignment \
    REFERENCE_SEQUENCE=../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta \
    UNMAPPED_BAM=$read_bam \
    ALIGNED_BAM="${out_loc}_aligned.sorted.bam" \
    OUTPUT="${out_loc}_merged.bam" \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false \
    QUIET=true

#newer function: tags overlaps with exons and introns - better mapping of exon-exon junctions
"${DropSeq_path}/TagReadWithGeneFunction" \
    I="${out_loc}_merged.bam" \
    O="${out_loc}_star_gene_exon_tagged.bam" \
    ANNOTATIONS_FILE=../meta_data/GCF_000001635.26_GRCm38.p6_genomic.gtf \
    2>/dev/null

mkdir "${align_out_folder}/bead_errors"

#detect bead synthesis errors
"${DropSeq_path}/DetectBeadSubstitutionErrors" \
    I="${out_loc}_star_gene_exon_tagged.bam" \
    O="${align_out_folder}/bead_errors/${group_num}_my_clean_subtitution.bam" \
    OUTPUT_REPORT="${align_out_folder}/bead_errors/${group_num}_my_clean.substitution_report.txt" \
    QUIET=true

#detect bead synthesis errors: need to add primer sequences to this
"${DropSeq_path}/DetectBeadSynthesisErrors" \
    I="${align_out_folder}/bead_errors/${group_num}_my_clean_subtitution.bam" \
    O="${align_out_folder}/bead_errors/${group_num}_my_clean.bam" \
    REPORT="${align_out_folder}/bead_errors/${group_num}_my_clean.indel_report.txt" \
    OUTPUT_STATS="${align_out_folder}/bead_errors/${group_num}_my.synthesis_stats.txt" \
    SUMMARY="${align_out_folder}/bead_errors/${group_num}_my.synthesis_stats.summary.txt" \
    QUIET=true

"${DropSeq_path}/BamTagHistogram" \
    I="${out_loc}_star_gene_exon_tagged.bam" \
    O="${out_loc}_out_cell_readcounts.txt.gz" \
    TAG=XC \
    QUIET=true

#use DigitalExpression to get gene expression matrix
"${DropSeq_path}/DigitalExpression" \
    I="${out_loc}_star_gene_exon_tagged.bam" \
    O="${out_loc}_out_gene_exon_tagged.dge.txt.gz" \
    SUMMARY="${out_loc}_out_gene_exon_tagged.dge.summary.txt" \
    NUM_CORE_BARCODES=100 \
    MIN_NUM_GENES_PER_CELL=$min_genes \
    QUIET=true

gunzip "${out_loc}_out_gene_exon_tagged.dge.txt.gz"

echo -e "${path_to_DGE}" | awk '{print $1 "\t" $2 "\t" $3}' >> "../${output_dir}/DGE-paths.txt"

done < "$input_file"
