import os 
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description="test_wrapper.py") #specify the name of the script
    parser.add_argument("-m", "--metadata", #add metadata file argument
    help="metadata with paired fastq names, sample group, number, and output location",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments 
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
sample_meta = arguments.metadata #load fastq directory

#sample_meta - location of sample txt
drop_seq_qc_command = "./drop_seq_qc.sh -i " + sample_meta

os.system(drop_seq_qc_command)

#os.system("wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz -P ../sample_data")
#os.system("wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz -P ../sample_data")

#os.system("gunzip ../sample_data/GCF_000001635.26_GRCm38.p6_genomic.fna.gz")
#os.system("gunzip ../sample_data/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz")

star_quant_command = "./drop_seq_mm10.sh -i " + sample_meta + " -o align_quant_output_test_mm10"

os.system(star_quant_command)

#seurat steps to do within the wrapper

