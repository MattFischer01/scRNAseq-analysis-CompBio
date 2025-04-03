#add argparse and call this from drop_seq_p2.sh
#for now: not performing it bc we already have a pre-made index

from Bio import SeqIO

fasta_file = "/home/project6/dropseq-qc-output/sample_data_qc/naive_unaligned_mc_tagged_polyA_filtered.fastq"  # Replace with your FASTA file

lengths = list()

for record in SeqIO.parse(fasta_file, "fastq"):
    lengths.append(len(record.seq))

max_read = max(lengths)

file = open("/home/project6/align_quant_output/max_read.txt", "w")
file.write(str(max_read - 1))
file.close()
