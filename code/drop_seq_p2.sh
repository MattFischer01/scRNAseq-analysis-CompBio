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

while IFS=$'\t' read -r read1 read2 output group number; do
    qc_output="${output}/${group}/${group}_${number}"
    align_output="${output_dir}/${group}/${group}_${number}"
    log_output="${output_dir}/${group}/${group}_${number}/log"
    log_file="${output_dir}/${group}_${number}_log.txt"

    mkdir -p "$align_output" "$log_output"

    echo "Processing: $read1 and $read2 > QC output: $align_output, Logs: $log_output"

#create a genome index for STAR
#we need to start a new directory to hold the indices

overhang = $(cat max_read.txt)

if ! test -d "genome_index"; then
  mkdir genome_index

  /home/project6/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
  --runThreadN 2 \
  --runMode genomeGenerate \
  --genomeDir genome_index \
  --genomeFastaFiles /home/project6/metadata/GCF_000001635.26_GRCm38.p6_genomic.fna \
  --sjdbGTFfile /home/project6/metadata/GCF_000001635.26_GRCm38.p6_genomic.gtf \
  --sjdbOverhang $overhang
fi

cd $output_dir

#before this: calculate maximum read length in fastq file using read_length.py - in this case, max is 80, so overhang is 79
#/home/project6/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
#--runThreadN 2 \
#--runMode genomeGenerate \
#--genomeDir /home/project6/genome_index_new \
#--genomeFastaFiles /home/project6/dropseq_metadata/metadata_mouse/GCF_000001635.26_GRCm38.p6_genomic.fna \
#--sjdbGTFfile /home/project6/dropseq_metadata/metadata_mouse/GCF_000001635.26_GRCm38.p6_genomic.gtf \
#--sjdbOverhang 79

#let's perform an alignment using STAR
/home/project6/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
--genomeDir genome_index \
--readFilesIn "$qc_output/${group}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.fastq" \
--outFileNamePrefix star

#sort output of alignment in queryname order
java -Xmx4g -jar /home/project6/tools/picard.jar SortSam \
    I="${group}/${group}_${number}_starAligned.out.sam" \
    O="${group}/${group}_${number}_aligned.sorted.bam" \
    SO=queryname

#merge alignment BAM with unaligned 
java -Xmx4g -jar /home/project6/tools/picard.jar MergeBamAlignment \
    REFERENCE_SEQUENCE=/home/project6/metadata/GCF_000001635.26_GRCm38.p6_genomic.fna \
    UNMAPPED_BAM="$qc_output/${group}/${group}_${number}_unaligned_mc_tagged_polyA_filtered.bam" \
    ALIGNED_BAM="${group}/${group}_${number}_aligned.sorted.bam" \
    OUTPUT="${group}/${group}_${number}_merged.bam" \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false

#newer function: tags overlaps with exons and introns - better mapping of exon-exon junctions
/home/project6/tools/dropseq-3.0.2/TagReadWithGeneFunction \
    I="${group}/${group}_${number}_merged.bam" \
    O="${group}/${group}_${number}_star_gene_exon_tagged.bam" \
    ANNOTATIONS_FILE=/home/project6/metadata/GCF_000001635.26_GRCm38.p6_genomic.gtf

mkdir bead_errors

#detect bead synthesis errors
/home/project6/tools/dropseq-3.0.2/DetectBeadSubstitutionErrors \
    I="${group}/${group}_${number}_star_gene_exon_tagged.bam" \
    O="bead_errors/${group}_${number}_my_clean_subtitution.bam" \
    OUTPUT_REPORT="bead_errors/${group}_${number}_my_clean.substitution_report.txt"

#detect bead synthesis errors: need to add primer sequences to this
/home/project6/tools/dropseq-3.0.2/DetectBeadSynthesisErrors \
    I="bead_errors/${group}_${number}_my_clean_subtitution.bam" \
    O="bead_errors/${group}_${number}_my_clean.bam" \
    REPORT="bead_errors/${group}_${number}_my_clean.indel_report.txt" \
    OUTPUT_STATS="bead_errors/${group}_${number}_my.synthesis_stats.txt" \
    SUMMARY="bead_errors/${group}_${number}_my.synthesis_stats.summary.txt" \
    PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

/home/project6/tools/dropseq-3.0.2/BamTagHistogram \
    I="${group}/${group}_${number}_star_gene_exon_tagged.bam" \
    O="${group}/${group}_${number}_out_cell_readcounts.txt.gz" \
    TAG=XC

gunzip "${group}/${group}_${number}_out_cell_readcounts.txt.gz"

python3 /home/2025/estagaman/scRNAseq-analysis-CompBio/code/cell_barcodes.py -i "${group}/${group}_${number}_out_cell_readcounts.txt.gz" -o  "${group}/${group}_${number}"

#use DigitalExpression to get gene expression matrix
/home/project6/tools/dropseq-3.0.2/DigitalExpression \
    I="${group}/${group}_${number}_star_gene_exon_tagged.bam" \
    O="${group}/${group}_${number}_out_gene_exon_tagged.dge.txt.gz" \
    SUMMARY="${group}/${group}_${number}_out_gene_exon_tagged.dge.summary.txt" \
    NUM_CORE_BARCODES=100 \
    MIN_NUM_GENES_PER_CELL=500 \
    CELL_BC_FILE="${group}/${group}_${number}_barcodes.txt"

gunzip "${group}/${group}_${number}_out_gene_exon_tagged.dge.txt.gz"
