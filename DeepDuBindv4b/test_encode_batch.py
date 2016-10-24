import DreamDataInput #Pull data batches

import numpy as np
import os
import sys
import time


def main():
    '''
    Purpose: Save an entire epoch of training data in the HDF5 format to speed up training.
    This should increase the speed of training by eliminating all pysam/bigwig calls
    and data preprocessing compute costs.
    '''

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





    print data_input.num_train_examples

    output_file = 'test_batch.h5'

    #for _ in range(data_input.num_train_examples):
    for _ in range (1):
        #t0 = time.clock()
        (dna_seq_batch,
        dna_shape_batch,
        dnase_seq_batch,
        chip_labels) = data_input.train.pull_batch_train(1)

        #composite_batch = np.concatenate((dna_seq_batch,
        #                                  dna_shape_batch,
        #                                  dnase_seq_batch,
        #                                  chip_labels),axis =0)



        #print "Time to pull (s):",time.clock()-t0

        print "Sequence:",dna_seq_batch.shape
        print "Shape:",dna_shape_batch.shape
        print "Dnase I accessibility", dnase_seq_batch.shape
        print "Labels", chip_labels.shape
        #print "\n"



        data_mat = None
        #with h5py.File(output_file,'a') as hf:
        #    try:
        #        hf.create_dataset('data',data=data_mat,chunks =True)
        #    except:
        #        print "Error saving data"

        print "Squeeze shape",np.squeeze(dna_seq_batch).shape
        print "Flat shape", np.ravel(dna_seq_batch).shape
        print "Flat squeeze shape", np.ravel(np.squeeze(dna_seq_batch)).shape

        flat1d = np.ravel(np.squeeze(dna_seq_batch))
        restored = np.reshape(flat1d,(1,1,600,4))
        print "Restored shape:", restored.shape

        print "Restored values", restored[0,0,0:20,3]
        print "Orig values    ", dna_seq_batch[0,0,0:20,3] 

        shape_list=determine_shapes(dna_seq_batch,
                                    dna_shape_batch,
                                    dnase_seq_batch,
                                    chip_labels)

        data_1d = encode_nd_to_1d(dna_seq_batch,
                                    dna_shape_batch,
                                    dnase_seq_batch,
                                    chip_labels)


        rec_dna,rec_shape,rec_seq, rec_chip = decode_1d_to_nd(data_1d,shape_list)
        print "Recovered shape", rec_chip.shape,rec_chip

        

        
def determine_shapes(*argv):
    #Determine the shapes of numpy arrays
    shape_list = []
    for arg in argv:
        shape_list.append(arg.shape)
    return shape_list
    
def encode_nd_to_1d(dna_seq,dna_shape,dnase_seq,chip_labels):
    #Flattens a batch to 1-D for saving into an HDF5 file
    dna_seq_1d = np.ravel(np.squeeze(dna_seq))
    dna_shape_1d = np.ravel(np.squeeze(dna_shape))
    dnase_seq_1d = np.ravel(np.squeeze(dnase_seq))
    chip_labels_1d = np.ravel(np.squeeze(chip_labels))
    return np.concatenate((dna_seq_1d,dna_shape_1d,dnase_seq_1d,chip_labels_1d),axis=0)
    


def decode_1d_to_nd(input_1d, shape_list):
    #Decode a 1D matrix into the shapes in shape_list
    return_list = []
    start_ind = 0
     #Slice out return data in their proper shapes
    for shape in shape_list:
        elem_len= np.prod(shape)
        return_list.append(np.reshape(input_1d[start_ind:elem_len+start_ind],shape))
        start_ind += elem_len
    return return_list
        



if __name__=="__main__":
    main()
    
