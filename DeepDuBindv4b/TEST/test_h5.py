import numpy as np
import h5py
import os
import time

#test_fname = "fasta_example_one.fa.Roll.h5"
#test_fname = "/media/storage3-ntfs/ENCODE-DREAM/annotations/chromosomes/hdf5_v1/chr1.fa.MGW.h5"

data_name = os.path.splitext(os.path.basename(test_fname))[0]
print "Keyname = ", data_name
fopen_start = time.clock()
with h5py.File(test_fname,'r') as hf:
        
    start = time.clock()
    rec =  hf[data_name][400000:401000]
    rec =  hf[data_name][10:20]
    elapsed = time.clock()-start
    print rec
    print "Time taken:",elapsed
    print "Length of record:", len(rec)

print "Fopenclose:",time.clock()-fopen_start
