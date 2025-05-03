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
9. Other R packages needed: dplyr, patchwork, ggplot2, enrichR, HGNChelper, openxlsx, SingleR, and celldex, which can be downloaded in R using ```install.packages("[package name]")```

Streamline R Package Installation:
```
install.packages(c("dplyr", "MAST", "Seurat", "patchwork", "celldex", "SingleR", "enrichR", "HGNChelper","openxlsx", "ggplot2"))
```


# Pipeline in Detail:

## Step 0: Clone the repo, move into code, and download dependencies

`git clone https://github.com/MattFischer01/scRNAseq-analysis-CompBio.git`

`cd scRNAseq-analysis-CompBio/code`

`bash download_dependences.sh`

download_dependences.sh will download Drop-seq, Picard, and STAR dependencies into a tools directory that will be used in the wrapper. This code will also download the mouse metadata for STAR alignment. 


## OPTIONAL Test dataset: 

Once you are in the code directory and you have run the download_dependenes.sh, run this command: 

`python3 wrapper_script.py -m sampledata-to-qc.txt -t ../tools`

The -m parameter allows you to specify the path to the sample data containing paired-end fastq files and the associated condition and sample group

## Step 1: Fetch FASTQ Files

If there are FASTQ files you need to download from NCBI's SRA database, the accession IDs (SRR) must be in a text file (.txt) with each ID on a new line. Please see srr_ids.txt for reference. This file needs to be in the main repository directory as well (above the code directory). 

This step fetches the SRA files using fasterq-dump and then gzips them to save space. The files will be saved in a directory titled scRNA_SRA in the main repository directory. This step may take awhile depending on number of FASTQ files downloaded and size of each file. 

This is to be run prior to the main wrapper if you do not have paired FASTQ files already. If you want to see how this wrapper works, you can run it with the current SRR IDs inside the text file with the steps below.

```
# cd into the code directory
cd code/

# give execution permissions to fetch_fastq.sh
chmod +x fetch_fastq.sh

# run script with srr_ids.txt as argument
./fetch_fastq.sh ../srr_ids.txt
```
This will create a scRNA__SRA directory in the main directory where the srr_ids.txt is. You may see temp files being created in the code directory, but this does not mean the FASTQ files are downloading in the code directory, the end file will be located in the right directory and the temp file will be deleted from the code directory.

## Step 2: Running the wrapper.py script

Here is a generalized version of the command:
```
python3 wrapper_script.py -m <<your metadata file>> -t <<tools directory>> -i <<genome index directory>> -f <<minimum number of genes per cell>> 
```
If running our test data, use: 
```
python3 wrapper_script.py -m sampledata-to-qc.txt -t ../tools -i /home/project6/genome_index_new -f 50 
```

Our wrapper script includes the following parameters: 

**--m** = path to metadata file  
&nbsp;The metadata file should include a tab-separated line for each sample, formatted similarly to `sampledata-to-qc.txt`. Each line should contain:  
		&nbsp;&nbsp;1. Path to the forward FASTQ  
		&nbsp;&nbsp;2. Path to the reverse FASTQ  
		&nbsp;&nbsp;3. Path to the output directory for QC results  
		&nbsp;&nbsp;4. Group identifier  
		&nbsp;&nbsp;5. Sample number  
&nbsp;These values must be listed in this exact order.

**--t** = path to tools directory  
	&nbsp;If `download_dependencies.sh` was used, the tools will be available in `../tools`. Otherwise, specify the location where you downloaded them.

**--i** = path to genome index generated by STAR  
	&nbsp;If `"auto"` is specified, STAR will generate the genome index for you (this can take 1–2 hours).  
	&nbsp;To use a pre-generated index: Use `/home/project6/genome_index_new`  

**--f** = minimum number of genes required to include a cell in the final gene expression matrix  
	&nbsp;In the original paper we aimed to replicate, `f` was set to `500`. However, since the test data is much smaller, we recommend using an `f` of `50`.

The wrapper script will run through the following steps: 

### Step 2a: DropSeq - QC
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


#### Output:
Each read pair following qc will be outputted as 1 FASTQ (since with appended the cell and molecular barcode to read 2 and deleted read 1) with the following file name: ${condition}_${sample number}_unaligned_mc_tagged_polyA_filtered.fastq, which will continue in the pipeline to STAR alignment. 

### Step 2b: STAR - alignment 

a. **create a genome index**
   - depending on your organism, you will fetch a reference genome (.fasta) and genomic features file (.gtf) from NCBI for that organism
   - this step alone can take 1-2 hours
         - in this case, we pre-downloaded reference files for the mouse genome, and pre-generated an index for faster processing
b. **align fastq from QC step to your genome index**

### Step 2c: DropSeq - Quantification

a. **Sort alignment output**
   - ensures alignments from the same read stay together in the BAM file

b. **Merge BAM files**
	- Match the aligned BAM output from STAR to the barcodes in the original BAM file 
	- Ensures all our reads are still tagged appropriately

c. **Tag genes in each read**
	- Uses a .gtf file of genomic features (provided) to tag any reads that overlap with a gene exon or intron
	- Reads are tagged with the gene name, gene strand, and gene function
	
d. **Look for errors in cell barcode sequences**
	- Sometimes, substitutions occur in barcode sequences by chance, which makes it hard to trace reads back to the right cells
	- This step in the processes combines barcodes that appear due to error with results from the intended barcode

e. **Calculate the number of reads per cell**
	- Outputs a file out_cell_readcounts.txt, which serves as input for the next step
	
f. **Create Digital Expression Matrix**
	- Settings: min_num_genes_per_cell = 500, num_core_barcodes=100
	- Generates a gene expression count matrix, which can serve as input for Seurat downstream analysis


#### Output:
The code will create an output directory with the name you specify. In our wrapper script, this directory is called align_quant_output

Each FASTQ from the QC step will output one gene expression matrix with the following file name: ${condition}_${sample number}_out_gene_exon_tagged.dge.txt, which can be used as input for Seurat and MAST


## Step 3 Seurat Wrapper and Tutorial

To process the gene expression matrices and integrate all of the matrices into one Seurat object, the Seurat wrapper must be ran. To do this, the user needs to provide a text file titled "DGE-paths.txt" that has the condition, condition number (if you have many of the same conditions, number them), and the path to that gene expression matrix all on the same line, separated by a space. You can see an example in the repository. This text file needs to be in the same directory as the gene expression matrices.

To run the R script, the user needs to provide an argument for the path where the gene expression matrices are.

For example, this runs the script in the background, which is recommended as these steps can take time depending on how many and how large the files are:
```
nohup Rscript seurat_test.R [path/to/gene/expression/matrices] &
```

After running this, you should be able to then load the RData file that is in the main repository directory into your RStudio session:
```
load("seurat_post_int.RData")
```

The steps after this will be detailed in the scRNAseq tutorial: https://tjensen15.github.io/scRNAseq-analysis-CompBio/. The tutorial goes through applying PCA, clustering, differential gene expression analysis, cell type annotation, and pathway analysis. The database for ScType cell annotation is included in the main directory of the repository and is titled the same as in the tutorial.

### Running sample data

Sample data was not provided to run the integration wrapper due to Seurat cell size and feature constraints, but you can run through the tutorial. For running the sample data, the RData file includes the mice.combined object after integration. You should be able to follow most of the steps at the beginning of the Seurat tutorial. Pathway enrichment analysis may not be able to be run since the sample size is so small, and you may get errors in this step. 

In your R Session, load the sample data:
```
load("seurat_small_subset.RData")
```
**For Dr. Wheeler:** If you want to run the larger (still subsetted) Seurat object, here is the path to the RData file: home/2025/hjensen2/archive/scRNAseq-analysis-CompBio/GEO/test_seurat.RData
