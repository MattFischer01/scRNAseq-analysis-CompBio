from Bio import SeqIO

fasta_file = "/home/project6/dropseq-qc-output/naive_unaligned_mc_tagged_polyA_filtered.fastq"  # Replace with your FASTA file

lengths = list()

for record in SeqIO.parse(fasta_file, "fastq"):
    lengths.append(len(record.seq))

max_read = max(lengths)
print(max_read)
