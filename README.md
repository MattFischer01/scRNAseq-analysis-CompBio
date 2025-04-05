# scRNAseq-analysis-CompBio
Pipeline for scRNAseq analysis 

The main goal of this project is to reproduce the pipeline for scRNAseq analysis seen in figure 1 of “MAFG-driven astrocytes promote CNS inflammation”. Thus, we aim to provide a user-friendly step-by-step guide for processing and assembling single-cell paired-end reads, using Seurat to perform cluster analysis of cell populations and MAST for differentially expressed gene analysis, and visualizing the data to interpret results. Analyzing high-throughput scRNAseq data is highly dependent on the researcher and the questions being asked, so our long-term goal is to create a pipeline to process raw data into a manageable and user-friendly analysis tool for the researcher to ask their own questions. Please look at our wiki for more details about the project.

## Workflow Overview:
### Part 1:
1. Download Fastq RNA-seq files from SRA
2. Perform QC and prep read for alignment using Drop-seq
3. Use STAR to map reads to the mouse reference genome
4. Quantify Gene Expression using Drop-Seq

### Part 2:
5. Cluster cells with Seurat
6. Compare cell types'/clusters' DEGs using MAST
7. Pathway Analysis comparing clusters/cell types
8. Subset cells to those with astrocyte markers, perform 5-7 again. 

![scRNA analysis flowchart-2](https://github.com/user-attachments/assets/49613a04-b9c6-4a4c-97fb-ac970d8d3c5c)


## Dependencies and Data
The following tools must be downloaded and installed to operate the entire pipeline. 
1. SRA Toolkit: https://anaconda.org/bioconda/sra-tools. SRA Toolkit also requires that conda and bioconda be installed before downloading the SRA Toolkit with bioconda.
2. Conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html 
3. Bioconda: https://bioconda.github.io/ 
4. Drop-Seq: https://github.com/broadinstitute/Drop-seq 
5. Picard: https://github.com/broadinstitute/picard 
6. STAR: https://github.com/alexdobin/STAR 
7. Seurat - R Package: https://cran.r-project.org/web/packages/Seurat/index.html. This can be downloaded in R using ```install.packages("Seurat")```
8. MAST: https://github.com/RGLab/MAST. This can also be downloaded as an R package using ```install.packages("MAST")```
9. Other R packages needed: dplyr, patchwork, SingleR, and celldex, which can be downloaded in R using ```install.packages("[package name]")```

Streamline R Package Installation:
```
install.packages(c("dplyr", "MAST", "Seurat", "patchwork", "celldex", "SingleR"))
```


# Pipeline in Detail:

## Step 0: Clone the repo and move into it

`git clone https://github.com/MattFischer01/scRNAseq-analysis-CompBio.git`

`cd scRNAseq-analysis-CompBio`


## Step 1: Fetch FASTQ Files

If there are FASTQ files you need to download from NCBI's SRA database, the accession IDs (SRR) must be in a text file (.txt) with each ID on a new line. Please see srr_ids.txt for reference.

This step fetches the SRA files using fasterq-dump and then gzips them to save space. The files will be saved in a directory titled scRNA_SRA. This step may take awhile depending on number of FASTQ files downloaded and size of each file. 


## Step 2: DropSeq - QC
Our DropSeq QC pipeline is derived from the DropSeq Alignment Cookbook from the McCarroll Lab, the founders of DropSeq. You are welcome to explore the McCarroll lab's pdf in our github for further details outside this brief explanation. 

The QC step reads in a tab-separated txt file with the following columns:

read1 path | read2 path | output path | condition | sample number 
--- | --- | --- | --- | --- 
/path/to/SRR8878106_1.fastq.gz | /path/to/scRNA_SRA/naive/SRR8878106_2.fastq.gz | /path/to/dropseq-qc-output | naive | 01

Overall, the pipeline reads in the input txt file (reads-to-qc.txt) with that FASTQ paired-end reads obtained from Drop-seq sequencies libraries, where read 1 contains both the cell and molecular barcode that help decide which read is coming from which cell and read 2 contains the sequenced cDNA from the cell. The raw reads are converted to a Picard-qeryname-sorted BAM file using Picard FastqToSam program. The BAM file undergoes the following qc: 

a. **Tag cell barcodes**  
   - Extract the first 12 bases of Read 1.  
   - Copy this barcode to the genome read using tag 'XC'.

b. **Tag molecular barcodes**  
   - Extract the last 8 bases of Read 1.  
   - Copy this barcode to the genome read using tag 'XM'.  
   - Delete Read 1 after processing.

c. **Filter low-quality reads**  
   - Remove reads that do not meet the quality threshold.

d. **Trim 5' primer sequence**  
   - Perform SMART adapter trimming.

e. **Trim 3' polyA sequence**  
   - Remove polyA tails to reduce bias in downstream analysis.

f. **Convert BAM back to FASTQ for alignment**  
   - Extract sequences from BAM and convert them to FASTQ format for alignment.

### Example command:
`/path/to/drop_seq_qc.sh -i /path/to/reads-to-qc.txt`

### Output:
Each read pair following qc will be outputted as 1 FASTQ (since with appended the cell and molecular barcode to read 2 and deleted read 1) with the following file name: ${condition}_${sample number}_unaligned_mc_tagged_polyA_filtered.fastq, which will continue in the pipeline to STAR alignment. 