#download STAR to /home/project6/tools (tools folder)

#create a genome index for STAR
#we need to start a new directory to hold the indices

mkdir /home/project6/genome_index
cd /home/project6/align_quant_output/

#before this: calculate maximum read length in fastq file using read_length.py - in this case, max is 80, so overhang is 79
/home/project6/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
--runThreadN 2 \
--runMode genomeGenerate \
--genomeDir /home/project6/genome_index \
--genomeFastaFiles /home/project6/dropseq_metadata/mm10/mm10.fasta \
--sjdbGTFfile /home/project6/dropseq_metadata/mm10/mm10.gtf \
--sjdbOverhang 79

#let's perform an alignment using STAR
/home/project6/STAR/STAR-2.7.11b/bin/Linux_x86_64/STAR \
--genomeDir /home/project6/genome_index \
--readFilesIn /home/project6/dropseq-qc-output/naive_unaligned_mc_tagged_polyA_filtered.fastq \
--outFileNamePrefix star

#sort output of alignment in queryname order
java -Xmx4g -jar /home/project6/tools/picard.jar SortSam \
    I=starAligned.out.sam \
    O=aligned.sorted.bam \
    SO=queryname

#merge alignment BAM with unaligned 
java -Xmx4g -jar /home/project6/tools/picard.jar MergeBamAlignment \
    REFERENCE_SEQUENCE=/home/project6/dropseq_metadata/mm10/mm10.fasta \
    UNMAPPED_BAM=/home/project6/dropseq-qc-output/unaligned_mc_tagged_polyA_filtered.bam \
    ALIGNED_BAM=aligned.sorted.bam \
    OUTPUT=merged.bam \
    INCLUDE_SECONDARY_ALIGNMENTS=false \
    PAIRED_RUN=false

#newer function: tags overlaps with exons and introns - better mapping of exon-exon junctions
/home/project6/tools/dropseq-3.0.2/TagReadWithGeneFunction \
    I=merged.bam \
    O=star_gene_exon_tagged.bam \
    ANNOTATIONS_FILE=/home/project6/dropseq_metadata/mm10/mm10.gtf

mkdir bead_errors

#detect bead synthesis errors
/home/project6/tools/dropseq-3.0.2/DetectBeadSubstitutionErrors \
    I=star_gene_exon_tagged.bam \
    O=bead_errors/my_clean_subtitution.bam \
    OUTPUT_REPORT=bead_errors/my_clean.substitution_report.txt

#detect bead synthesis errors: need to add primer sequences to this
/home/project6/tools/dropseq-3.0.2/DetectBeadSynthesisErrors \
    I=bead_errors/my_clean_subtitution.bam \
    O=bead_errors/my_clean.bam \
    REPORT=bead_errors/my_clean.indel_report.txt \
    OUTPUT_STATS=bead_errors/my.synthesis_stats.txt \
    SUMMARY=bead_errors/my.synthesis_stats.summary.txt \
    PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

#use DigitalExpression to get gene expression matrix
/home/project6/tools/dropseq-3.0.2/DigitalExpression \
    I=star_gene_exon_tagged.bam \
    O=out_gene_exon_tagged.dge.txt.gz \
    SUMMARY=out_gene_exon_tagged.dge.summary.txt \
    NUM_CORE_BARCODES=100 \
    MIN_NUM_GENES_PER_CELL=500
