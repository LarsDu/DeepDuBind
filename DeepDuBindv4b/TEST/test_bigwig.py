
import pyBigWig
import numpy as np
import os
import time

#https://github.com/dpryan79/pyBigWig#installation
#wiggles_dir = ("../ENCODE_DREAM_LINK"+os.sep+"DNASE"+os.sep+"fold_coverage_wiggles")

#test_fname = "DNASE.GM12878.fc.signal.bigwig"

#fname = os.path.join(wiggles_dir,test_fname)

#start =time.clock()
#bw = pyBigWig.open(fname)
#bw= pyBigWig.open(fname)
#print bw.header()
#print time.clock()-start

#print(bw.chroms("chr1"))
    #print bw.header()
    
#beg=135013200
#end=135018000
#beg = 2
#end = 20
#nbins = int((end-beg)/4)
#print nbins
#print bw.stats('chr2',beg,end,type="mean",nBins=nbins)
#bw.close()    

print '\n\n\n\n\n'

wiggles_dir = ("../ENCODE_DREAM_LINK"+os.sep+"DNAshape"+
               os.sep+"BigWig")

test_fname = 'hg19.MGW.wig.bw'


fname = os.path.join(wiggles_dir,test_fname)
bw = pyBigWig.open(fname)
print fname
#print bw.chroms()
#print bw.header()
print "Max",bw.header()['maxVal']
print "Min",bw.header()['minVal']
print '\n'


mean_val = bw.header()['sumData']/bw.header()['nBasesCovered']
print 'Mean value:',mean_val
contig='chr22'
start= 1230000
end = 1230100
t0 = time.clock()
A = bw.values(contig,start,end)
print "Time A:", time.clock()-t0
#print A
t0 = time.clock()
B = bw.stats(contig,start,end,type='mean',nBins = 100)
print "Time B:", time.clock()-t0

#A = np.nan_to_num(A)
#print A

bw.close()

            
    
