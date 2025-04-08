#example command: python wrapper_script.py -m sampledata-to-qc.txt -t /home/project6/tools_all
import os 
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description="wrapper_mm10.py") #specify the name of the script
    parser.add_argument("-m", "--metadata", #add metadata file argument
    help="metadata with paired fastq names, sample group, number, and output location",
    required=True)

    parser.add_argument("-t", "--tools", #add metadata file argument
    help="directory containing tools Picard, Drop-Seq, and STAR",
    required=True)

    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
sample_meta = arguments.metadata #load file with paired-end reads, sample group, and output location
tools_dir = arguments.tools #load tools directory

#sample_meta - location of sample txt
drop_seq_qc_command = "./drop_seq_qc.sh -i " + sample_meta + " -t " + tools_dir

os.system(drop_seq_qc_command)

os.system("wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz -P ../sample_data")
os.system("wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz -P ../sample_data")
os.system("mv ../sample_data/GCF_000001635.26_GRCm38.p6_genomic.fna.gz ../sample_dataGCF_000001635.26_GRCm38.p6_genomic.fasta.gz")

os.system("gunzip ../sample_data/GCF_000001635.26_GRCm38.p6_genomic.fasta.gz")
os.system("gunzip ../sample_data/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz")

star_quant_command = "./drop_seq_p2.sh -i " + sample_meta + " -o align_quant_output -t " + tools_dir

os.system(star_quant_command)

#seurat steps to do within the wrapper

