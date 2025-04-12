#!/bin/bash

#download STAR to /home/project6/tools (tools folder)
#Example command to run: nohup /home/2025/mfischer10/scRNAseq-analysis-CompBio/code/drop_seq_p2.sh -i /home/2025/mfischer10/scRNAseq-analysis-CompBio/code/reads-to-qc.txt  -o align-quant-output -t /home/project6/tools_all &

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

if [[ "$5" == "-t" || "$3" == "--tools" ]]; then #take in tools directory as tools_dir
    tools_dir="$6"
else
    echo "Unknown parameter passed: $6"
    exit 1
fi

STAR=$(find ${tools_dir} -maxdepth 1 -name "*STAR*" -exec basename {} \;) #search in tools_dir for name containing STAR
DropSeq=$(find ${tools_dir} -maxdepth 1 -name "*dropseq*" -exec basename {} \;) #search in tools_dir for name containing dropseq
Picard=$(find ${tools_dir} -maxdepth 1 -name "*picard*" -exec basename {} \;) #search in tools_dir for name containing picard

STAR_path="${tools_dir}/${STAR}/bin/Linux_x86_64/STAR" #construct full STAR path
DropSeq_path="${tools_dir}/${DropSeq}" #construct full DropSeq path
Picard_path="${tools_dir}/${Picard}" #construct full Picard path

while IFS=$'\t' read -r read1 read2 output group number; do
    qc_output="${output}/${group}/${group}_${number}" #location to pull QC results from
    align_out_folder="../${output_dir}/${group}/${group}_${number}" #folder to put Align-Quant results
    log_output="../${output_dir}/${group}/${group}_${number}/log" #log directory
    log_file="../${output_dir}/${group}_${number}_log.txt" #log file name
    read_in="${qc_output}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.fastq" #file from QC to read into STAR
    read_bam="${qc_output}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.bam" #file from QC to read into SortSam
    out_loc="${align_out_folder}/${group}_${number}" #output location including prefix
    group_num="${group}_${number}" #just prefix

    mkdir -p "${align_out_folder}" "$log_output" #make the align-quant output directory and log directory

    echo "Processing: $read1 and $read2 > alignment output: $align_output, Logs: $log_output"

#####This code generates a genome index. It takes 1-2 hours, so I have commented it out for now########

#overhang = $(cat max_read.txt)

#if ! test -d "../meta_data/genome_index"; then
  #mkdir ../meta_data/genome_index

  #$STAR_path \
  #--runThreadN 2 \
  #--runMode genomeGenerate \
  #--genomeDir ../meta_data/genome_index \
  #--genomeFastaFiles ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta \
  #--sjdbGTFfile ../meta_data/GCF_000001635.26_GRCm38.p6_genomic.gtf \
  #--sjdbOverhang 79
#fi

#let's perform an alignment using STAR: 
    #note: index is hard-coded for now, but the user should either supply the location to this as input or have the code generate it
$STAR_path \
--genomeDir /home/project6/genome_index_new \
--readFilesIn $read_in \
--outFileNamePrefix "${out_loc}_star"

#sort output of alignment in queryname order
java -Xmx4g -jar $Picard_path SortSam \
    I="${out_loc}_starAligned.out.sam" \
    O="${out_loc}_aligned.sorted.bam" \
    SO=queryname

#####creating metadata from the mouse reference files which can be used later
java -jar $Picard_path CreateSequenceDictionary \
    R=../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta \
    O=../meta_data/GCF_000001635.26_GRCm38.p6_genomic.dict \

#merge alignment BAM with unaligned
java -Xmx4g -jar $Picard_path MergeBamAlignment \
    REFERENCE_SEQUENCE=../meta_data/GCF_000001635.26_GRCm38.p6_genomic.fasta \
    UNMAPPED_BAM=$read_bam \
    ALIGNED_BAM="${out_loc}_aligned.sorted.bam" \
    OUTPUT="${out_loc}_merged.bam" \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false

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