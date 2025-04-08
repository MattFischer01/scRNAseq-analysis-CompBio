#example command: python test_wrapper.py -m sampledata-to-qc.txt -t /home/project6/tools_all

import os 
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description="test_wrapper.py") #specify the name of the script
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

os.system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz -P ../sample_data")
os.system("tar -xvzf ../sample_data/GSE63472_mm10_reference_metadata.tar.gz -C ../sample_data")

star_quant_command = "./drop_seq_mm10.sh -i " + sample_meta + " -o align_quant_output -t " + tools_dir

os.system(star_quant_command)

#seurat steps to do within the wrapper

