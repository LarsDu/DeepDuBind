import os
import gzip
import time

fname ='ENCODE_DREAM_LINK/annotations/train_regions.blacklistfiltered.bed.gz'
#ladder_fname=''
#train_fname=''

if os.path.splitext(fname)[1] == '.gz':
    print "confirmed"
with gzip.open(fname,'r') as f:
    lines = f.readlines()
    print lines[0].strip().split('\t')


    
def visit_file(fname):
    with gzip.open(fname,'r') as f:
        for line in f:
            contig, start, end = line.strip().split('\t')
            start = int(start)
            end = int(end)
    
