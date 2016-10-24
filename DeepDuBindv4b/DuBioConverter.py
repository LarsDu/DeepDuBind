import os
import glob
from subprocess import Popen,PIPE


def main():
    print "Testing"
    #test_glob('/home/ladu/Code/ENCODE_DREAM_LINK/ChIPseq/peaks/relaxed/test_dir')
#    narrowpeak_gz_to_bgz('/home/ladu/Code/ENCODE_DREAM_LINK/ChIPseq/peaks/relaxed/')
    narrowpeak_gz_to_bgz('/home/ladu/Code/ENCODE_DREAM_LINK/DNASE/peaks/relaxed/')

    #narrowpeak_gz_to_bigbed('/home/ladu/Code/ENCODE_DREAM_LINK/ChIPseq/peaks/relaxed/test_dir',
    #                        '../ENCODE_DREAM_LINK/annotations/hg19.chrom.sizes')

def test_glob(dir):
    os.chdir(dir)
    for fname in glob.glob('*narrowPeak.gz'):
        print os.path.basename(fname)
        print os.path.splitext(fname)[0]
    
def tsv_gz_to_bgz_dir(tsv_dir,sorted_fname_ext ='.sorted',do_remove_original =True):
    os.chdir(tsv_dir)
    for fname in glob.glob("*tsv.gz"):
        tsv_to_bgz(fname,sorted_fname_ext,do_remove_original)
        
def tsv_to_bgz(fname,sorted_fname_ext='.sorted',do_remove_original=False):

    print "Filename is",fname
    ext = os.path.splitext(fname)[1]
    bgz_fname = os.path.splitext(fname)+sorted_fname_ext+'.bgz'
    if os.path.isfile(bgz_fname):
        print 'File',bgz_fname,'already exists. Skipping conversion...'
    else:
        if ext == '.gz' or ext == '.gzip':
            print "Decompressing",fname
            Popen(['gzip','-df',fname],stdout=PIPE)
            unzip_fname = os.path.splitext(fname)[0]
        elif ext == '.tsv':
            unzip_fname = fname
        sorted_fname = unzip_fname+sorted_fname_ext
        print "Sorting",unzip_fname
        with open(sorted_fname,'w') as sf:
            Popen(['sort', '-k1,1', '-k2,2n', unzip_fname],stdout=sf)

        print "Bgzip compressing",sorted_fname
        with open(bgz_fname,'w') as bgzf:
            Popen([ 'bgzip', '-c', sorted_fname],stdout=bgzf)
        print "Tabix indexing",bgz_fname
        Popen(['tabix', '-p', 'bed', bgz_fname],stdout=PIPE)
        if do_remove_original==True:
            Popen(['rm', '-f', fname],stdout=PIPE)


def narrowpeak_gz_to_bgz_dir(npeaks_dir,sorted_fname_ext ='.sorted',do_remove_original =True):
    print npeaks_dir
    os.chdir(npeaks_dir)
    for fname in glob.glob("*narrowPeak.gz"):
        narrowpeak_to_bgz(fname,sorted_fname_ext,do_remove_original)

def narrowpeak_gz_bgz_dir(npeaks_dir,sorted_fname_ext ='.sorted',do_remove_original =True):
    print npeaks_dir
    os.chdir(npeaks_dir)
    for fname in glob.glob("*narrowPeak"):
        narrowpeak_to_bgz(fname,sorted_fname_ext,do_remove_original)

                
def narrowpeak_to_bgz(fname,sorted_fname_ext='.sorted',do_remove_original=False):
    print "Filename is",fname
    
    ext = os.path.splitext(fname)[1]
    bgz_fname = os.path.splitext(fname)[0]+sorted_fname_ext+'.bgz'
    if os.path.isfile(bgz_fname):
        print 'File',bgz_fname,'already exists. Skipping conversion...'
    else:
        if ext == '.gz' or ext == '.gzip':
            print "Decompressing",fname
            Popen(['gzip','-df',fname],stdout=PIPE)
            unzip_fname = os.path.splitext(fname)[0]
        elif ext == '.narrowPeak':
            unzip_fname = fname
        sorted_fname = unzip_fname+sorted_fname_ext
        print "Sorting",fname
        with open(sorted_fname,'w') as sf:
            Popen(['sort', '-k1,1', '-k2,2n', fname],stdout=sf)
        print "Bgzip compressing",sorted_fname
        with open(bgz_fname,'w') as bgzf:
            Popen([ 'bgzip', '-c', sorted_fname],stdout=bgzf)
        print "Tabix indexing",bgz_fname
        Popen(['tabix', '-p', 'bed', bgz_fname],stdout=PIPE)
        if do_remove_original==True:
            Popen(['rm', '-f', fname],stdout=PIPE)


def narrowpeak_gz_to_bigbed_dir(bb_dir,
                            chrom_sizes_loc='.',
                            bedToBigBed_loc = './bedToBigBed',
                            sorted_fname_ext ='.sorted',do_remove_original =False):

    #Note this appears to be broken for the contest data files due to lack of bed headers
    os.chdir(bb_dir)
    for fname in glob.glob("*narrowPeak.gz"):
        print "Filename is",fname
        unzip_fname = os.path.splitext(fname)[0]
        sorted_fname = unzip_fname+sorted_fname_ext
        bb_fname = sorted_fname+'.bb'
        if os.path.isfile(bb_fname):
            print 'File',bb_fname,'already exists. Skipping conversion...'
        else:
            print "Decompressing",fname
            Popen(['gzip','-df',fname],stdout=PIPE)
        print "Sorting",unzip_fname
        with open(sorted_fname) as sf:
            Popen(['sort', '-k1,1', '-k2,2n', unzip_fname],stdout=sf)
        print "bedToBigBed converting",sorted_fname
        proc = Popen([bedToBigBed_loc, sorted_fname,chrom_sizes_loc,bb_fname],stdout=PIPE)
        for line in proc.stdout:
            print line
        if do_remove_original==True:
            Popen(['rm', '-f', fname],stdout=PIPE)

                
                




    
if __name__ =="__main__":
    main()
