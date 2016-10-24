import time
import pysam

fname = "../ENCODE_DREAM_LINK/ChIPseq/peaks/relaxed/ChIPseq.GM12878.EGR1.relaxed.narrowPeak.gz" 
npf = pysam.TabixFile(fname)
t0 = time.clock()
records = npf.query('chr10,100000,100500')
for rec in records:
    print rec
