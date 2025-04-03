import argparse
#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="cell_barcodes.py") #specify the name of the script
    parser.add_argument("-o", "--output", #add metadata file argument
    help="output prefix",
    required=True)
    parser.add_argument("-i", "--input", #add directory with test data argument
    help="input file",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments 
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
infile = arguments.input #load metadata
outfile = arguments.output #load fastq directory

import pandas as pd

barcodes = pd.read_csv("out_cell_readcounts.txt", sep='\t', skiprows=1, header=None, names=['count', 'barcode'])

barcodes = barcodes.drop('count', axis=1)

out_loc = outfile + "_barcodes.txt"

barcodes.to_csv(out_loc, header=None, index=None, sep=' ', mode='a')
