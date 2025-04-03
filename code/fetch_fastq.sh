#!/bin/bash

# check if an argument (file name) is provided
if [ "$#" -ne 1 ]; then # if the number of arguments passed is not exactly 1
    echo "Usage: $0 <srr_ids_file>" # print usage instructions ($0 is the script name)
    exit 1 # exit script with error
fi

SRR_FILE="$1" # assign first argument to the variable SRR_FILE

mkdir -p SRR  # '-p' ensures no error if the directory already exists

# read SRR IDs from the provided file and process them sequentially to ensure each is processed
while read -r id; do
    echo "Processing $id..."
    
    # run fasterq-dump and output files directly into the SRR directory
    fasterq-dump "$id" -O SRR && gzip SRR/"$id"*.fastq  

    if [ $? -eq 0 ]; then # loop states whether or not the file is processing
        echo "Successfully processed $id"
    else
        echo "Error processing $id"
    fi
done < "$SRR_FILE"
