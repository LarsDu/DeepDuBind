import DreamDataInput #Pull data batches

import numpy as np
import os
import sys
import time




base_dir = '/home/ladu/Code/ENCODE_DREAM_LINK/'
train_coords_file=base_dir+'train_test/train_regions.blacklistfiltered.merged.bed'
test_coords_file=base_dir+'train_test/train_regions.blacklistfiltered.merged.bed'
eval_coords_file=base_dir+'train_test/train_regions.blacklistfiltered.merged.bed'

genome_file=base_dir+'hg19_genome/hg19.genome.fa'
genome_chromosome_sizes=base_dir+'annotations/hg19.chrom.sizes'
genome_annotation_file=base_dir+'annotations/gencode.v19.annotation.gff3'
dna_shape_dir=base_dir+'DNAshape/BigWig'
dnase_seq_dir=base_dir+'DNASE/fold_coverage_wiggles/'
chip_seq_file=base_dir+'ChIPseq/labels/EGR1.train.labels.tsv.bgz'
test_cell_type='K562'
eval_cell_type='K562'
max_train_steps=45000
learning_rate=1e-4
batch_sisze = 48


data_input = DreamDataInput.DreamDataLoader(
                                'train',
                                train_coords_file,
                                test_coords_file,
                                eval_coords_file,
                                None,
                                test_cell_type,
                                eval_cell_type,
                                genome_file,
                                genome_chromosome_sizes,
                                dna_shape_dir,
                                dnase_seq_dir,
                                None,
                                genome_annotation_file,
                                chip_seq_file,
                                dna_window = 600,
                                chrom_window = 600*6,
                                k_folds = 0,
                                k_validation_test_frac=.5
                                )




test_shape = data_input.dna_shape_data.load_single_coords_bw('chr10',150000,150100)

print test_shape

for _ in range(1):
    t0 = time.clock()
    (dna_seq_batch,
    dna_shape_batch,
    dnase_seq_batch,
    chip_labels) = data_input.train.pull_batch_train(1)
    print "Time to pull (s):",time.clock()-t0


    print "DNA Sequence Shape:",dna_seq_batch.shape
    print "DNA Shape Shape:",dna_shape_batch.shape
    print "Dnase I accessibility Shape", dnase_seq_batch.shape
    print "Labels Shape", chip_labels.shape
    print "DNA Sequence:",dna_seq_batch
    print "DNA Shape:",dna_shape_batch
    print "Dnase I accessibility", dnase_seq_batch
    print "Labels", chip_labels
    print "\n"
