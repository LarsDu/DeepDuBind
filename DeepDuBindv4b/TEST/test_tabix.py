import tabix
import dutabix
import time
from subprocess import Popen, PIPE
import sys
import pysam

#fname = "ENCODE_DREAM_LINK/ChIPseq/labels/TAF1.train.labels.tsv.bgz"

fname = 'ENCODE_DREAM_LINK/ChIPseq/labels/ARID3A.train.labels.tsv.bgz'
print "Timing pytabix:"

tb = tabix.open(fname)
t0 = time.clock()
records = tb.query('chr10',100000,100500)
print time.clock()-t0
for record in records:
    print record

print "Timing pysam"

tb2 = pysam.TabixFile(fname)
t0 = time.clock()
records = tb.query('chr10',100000,100500)
print time.clock()-t0
for record in records:
    print record.strip().split()






#for record in records:
#    print record

def tabix_query(filename, chrom, start, end):
        """Call tabix and generate an array of strings for each line it returns."""
        query = '{}:{}-{}'.format(chrom, start, end)
        process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
        for line in process.stdout:
            yield line.strip().split()
print "Tabix_query"        
t0=time.clock()
records = tabix_query(fname,'chr10',575000,575025)
print time.clock()-t0


for record in records:
    print record

print "Testing dutabix popen"
t0=time.clock()
records = dutabix.query_popen(fname,'chr10',575000,575100)
sys.stdout.flush()
for rec in records:
    print rec
print time.clock()-t0
#t0=time.clock()
#qt = dutabix.query_commands(fname,'chr10',575000,575200)
#print "Test", qt[0]
#print time.clock()-t0
#Conclusion: pytabix has some insanely slow overhead
