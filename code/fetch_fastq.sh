#!/bin/bash

# check if an argument (file name) is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <srr_ids_file>"
    exit 1
fi

# get the absolute path of the provided SRR text file
SRR_FILE="$(realpath "$1")"

# define the output directory
OUTPUT_DIR="$(dirname "$SRR_FILE")/scRNA_SRA"

# create the output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# read SRR IDs from the provided file and process them 
while read -r id; do
    echo "Processing $id..."
    
    # run fasterq-dump and output files directly into the scRNA_SRA directory
    fasterq-dump "$id" -O "$OUTPUT_DIR" && gzip "$OUTPUT_DIR/$id"*.fastq  

    # this if then else statements are to echo whether or not the srr id is being processed
    if [ $? -eq 0 ]; then
        echo "Successfully processed $id"
    else
        echo "Error processing $id"
    fi
done < "$SRR_FILE"
