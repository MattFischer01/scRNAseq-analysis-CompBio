#example command: python wrapper_script.py -m sampledata-to-qc.txt -t /home/project6/tools_all -i /home/project6/genome_index_new -f 50
import os
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description="wrapper_script.py") #specify the name of the script
    parser.add_argument("-m", "--metadata", #add metadata file argument
    help="metadata with paired fastq names, sample group, number, and output location",
    required=True)

    parser.add_argument("-t", "--tools", #add metadata file argument
    help="directory containing tools Picard, Drop-Seq, and STAR",
    required=True)

    parser.add_argument("-i", "--index", #add metadata file argument
    help="directory containing the index. If 'auto' is provided, STAR will generate the index",
    required=True)

    parser.add_argument("-f", "--filter_by_genes", #add metadata file argument
    help="directory containing the index. If 'auto' is provided, STAR will generate the index",
    required=True)

    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
sample_meta = arguments.metadata #load file with paired-end reads, sample group, and output location
tools_dir = arguments.tools #load tools directory
index_loc = arguments.index #load index directory - if index = "auto", STAR will generate the index
min_genes = arguments.filter_by_genes

#sample_meta - location of sample txt
drop_seq_qc_command = "./drop_seq_qc.sh -i " + sample_meta + " -t " + tools_dir

os.system(drop_seq_qc_command)

star_quant_command = "./drop_seq_p2.sh -i " + sample_meta + " -o align-quant-output -t " + tools_dir + " -x " + index_loc + " -g " + min_genes

os.system(star_quant_command)

#seurat steps to do within the wrapper