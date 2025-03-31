#Step 0: Extract metadata used in alignment:
#The Dropseq pre-made metadata is from an earlier Mouse Genome (GRCm38) https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.26/ 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz
tar -xvzf GSE63472_mm10_reference_metadata.tar.gz 

#Step0.1: The MAFG paper uses a different mouse reference genome (GRCm38.p6) https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index 
#I am pulling from this ftp site: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/ 
#genome 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
# then convert to fasta:
mv GCF_000001635.26_GRCm38.p6_genomic.fna.gz GCF_000001635.26_GRCm38.p6_genomic.fasta.gz
#gtf 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz


##Step1: Extract necessary tools:
wget https://github.com/broadinstitute/Drop-seq/releases/download/v3.0.2/dropseq-3.0.2.zip

unzip dropseq-3.0.2.zip 

wget https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar

#Creation of metadata: CreateSequenceDictionary 
java -jar /home/project6/tools/picard.jar CreateSequenceDictionary \
    R=/home/project6/dropseq_metadata/metadata_mouse/GCF_000001635.26_GRCm38.p6_genomic.fna.gz \
    O=/home/project6/dropseq_metadata/metadata_mouse/mouse.dict

/home/project6/tools/dropseq-3.0.2/ConvertToRefFlat \
    ANNOTATIONS_FILE=/home/project6/dropseq_metadata/metadata_mouse/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz\
    SEQUENCE_DICTIONARY=/home/project6/dropseq_metadata/metadata_mouse/mouse.dict \
    OUTPUT=/home/project6/dropseq_metadata/metadata_mouse/mouse.refFlat

/home/project6/tools/dropseq-3.0.2/ReduceGTF \
    SEQUENCE_DICTIONARY=/home/project6/dropseq_metadata/metadata_mouse/mouse.dict \
    GTF=/home/project6/dropseq_metadata/metadata_mouse/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz \
    OUTPUT=/home/project6/dropseq_metadata/metadata_mouse/mouse.reduced.gtf


#4:30pm, 03/29/24, this is taking far too long to figure out for the updated reference genome, let's continue with pre-made meta data. 