import subprocess
import commands
import sys

#This code was from https://github.com/slowkow/pytabix
#For some reason, these subprocess examples run vastly faster
#than the pytabix module.

def bgzip(filename):
    """Call bgzip to compress a file."""
    Popen(['bgzip', '-f', filename])

def index(filename,
    preset="gff", chrom=1, start=4, end=5, skip=0, comment="#"):
    """Call tabix to create an index for a bgzip-compressed file."""
    Popen(['tabix', '-p', preset, '-s', chrom, '-b', start, '-e', end,
                '-S', skip, '-c', comment])

def query_popen(filename, chrom, start, end):
        """Call tabix and generate an array of strings for each line it returns."""
        query = '{}:{}-{}'.format(chrom, start, end)
        process = subprocess.Popen(['tabix', filename, query],stdout=subprocess.PIPE)
        for line in process.stdout:
            yield line.strip().split()
            

def query_commands(filename,chrom,start,end):
    
    query = '{}:{}-{}'.format(chrom, start, end)
    output = commands.getoutput(' '.join(['tabix',filename,query]) )
    output = output.split('\n')
    for line in output:
        yield line.strip().split('\t')
    
