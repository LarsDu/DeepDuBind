import pysam
from Bio import SeqIO
import time
import pyfaidx

###pysam comparison
fasta_fname = "ENCODE_DREAM_LINK/annotations/hg19.genome.fa"
start = time.clock()
genome = pysam.FastaFile(fasta_fname)
print "samfetch load"
print str(time.clock()-start)


start=time.clock()
print genome.fetch("chr10",400000,400100)
print "samfetch"
print str(time.clock()-start)

###pyfaidx comparison
fasta_fname = "ENCODE_DREAM_LINK/annotations/hg19.genome.fa"
start = time.clock()
genome = pyfaidx.Fasta(fasta_fname)
print "pyfaidx load"
print str(time.clock()-start)


start=time.clock()
print genome['chr10'][400000:400100]
print "pyfaidx fetch"
print str(time.clock()-start)

###pyfaidx for getting DNAShapeR data
#fasta_fname = "ENCODE_DREAM_LINK/annotations/chromosomes/chr10.fa.HelT"
#fasta_fname = "fasta_example.fa.HelT.bytefix"
start = time.clock()
genome = pyfaidx.Fasta(fasta_fname)
print "pyfaidx raw extract"
print str(time.clock()-start)


start=time.clock()
print genome['chr10'][400000:400100]
print "pyfaidx raw extract"
print str(time.clock()-start)



###Biopython comparison
#Memory intensive non-indexed technique

#start = time.clock()
#genome_dict = SeqIO.index(fasta_fname,"fasta")
#print "biopython load"
#print str(time.clock()-start)


#start = time.clock()
#print genome_dict["chr10"].seq[400000:400100]
#print "biopython fetch"
#print str(time.clock()-start)



#Result: BioPython sequence fetch from fasta is slow as balls
# Loading and retrieving using pysam takes .00064 seconds whereas
# retrieving a sequence using biopython (with fasta completely loaded
# to memory takes 33 seconds to load and 2.3 seconds to retrieve)
