#example command: python wrapper_mm10.py -m sampledata-to-qc.txt -t /home/project6/tools_all

import os 
import argparse
import sys

def check_arg(args=None):
    parser = argparse.ArgumentParser(description="wrapper_mm10.py") #specify the name of the script
    parser.add_argument("-m", "--metadata", #add metadata file argument
    help="metadata with paired fastq names, sample group, number, and output location",
    required=True)

    parser.add_argument("-t", "--tools", #add tools directory argument
    help="directory containing tools Picard, Drop-Seq, and STAR",
    required=True)

    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
sample_meta = arguments.metadata #load file with paired-end reads, sample group, and output location
tools_dir = arguments.tools #load tools directory

#construct dropseq_qc command: -i is the sample metadata file, -t is the directory holding STAR, Picard, and DropSeq
drop_seq_qc_command = "./drop_seq_qc.sh -i " + sample_meta + " -t " + tools_dir
os.system(drop_seq_qc_command)

#Download mouse metadata (mm10) if it doesn't already exist in the sample_data folder
file_path = "../sample_data/GSE63472_mm10_reference_metadata.tar.gz"
mm10_path = "../sample_data/mm10/mm10.dict"
file_url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz"
to_unzip = "../sample_data/GSE63472_mm10_reference_metadata.tar.gz"

# Check if the file already exists
if os.path.isfile(file_path):
    print("mm10 reference file already exists.")
else:
    print("mm10 reference file does not exist. Proceeding with download...")
    
    os.system(f"wget {url} -P ../sample_data")
    print(f"File downloaded to: {file_path}")

#Check if it has already been unzipped
if os.path.isfile(file_path):
    print("mm10 reference file already unzipped.")
else:
    print("mm10 reference file needs to be unzipped. Unzipping...")
    
    os.system(f"tar -xvzf {to_unzip} -C ../sample_data")
    print(f"File downloaded to: {file_path}")


#construct STAR command: -i is the sample metadata file, -o is the name you would like the output folder to be called, -t is the tools directory
star_quant_command = "./drop_seq_mm10.sh -i " + sample_meta + " -o align_quant_output -t " + tools_dir
os.system(star_quant_command)

#NEXT: seurat steps to do within the wrapper

