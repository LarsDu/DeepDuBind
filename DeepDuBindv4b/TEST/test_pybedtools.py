import numpy as np
import pybedtools
from pybedtools import BedTool
import os
import DuBioTools as dbt

#https://www.biostars.org/p/158849/
base_dir = 'ENCODE_DREAM_LINK'+os.sep
coord_file = 'annotations'+os.sep+'test_regions.blacklistfiltered.bed.gz'
genome_file = 'annotations'+os.sep+'hg19.genome.fa'

coords = BedTool(base_dir+coord_file)
genome = BedTool(base_dir+genome_file)

extract = coords.sequence(fi=genome)

print (open(extract.seqfn).read())
