import pandas as pd

barcodes = pd.read_csv("out_cell_readcounts.txt", sep='\t', skiprows=1, header=None, names=['count', 'barcode'])

barcodes = barcodes.drop('count', axis=1)

barcodes.to_csv("barcodes.txt", header=None, index=None, sep=' ', mode='a')
