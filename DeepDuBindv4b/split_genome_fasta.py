#from Bio import SeqIO
#from Bio.Seq import Seq
import DuBioTools as dbt

#This script was written in tandem with split_genome_fasta.R to circumvent
# memory limitation on the DNAshapeR R package.
# This script should generate a fasta file for each chromosome of a given genome

#fname = "fasta_example.fa"
fname = './ENCODE_DREAM_LINK/annotations/hg19.genome.fa'

dbt.split_fasta(fname)
