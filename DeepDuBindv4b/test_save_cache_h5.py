import DreamDataInput

'''
Purpose: Save an entire epoch of training data in the HDF5 format to speed up training.
This should increase the speed of training by eliminating all pysam/bigwig calls
and data preprocessing compute costs.
'''

base_dir = '/home/ladu/Code/ENCODE_DREAM_LINK/'
train_coords_file=base_dir+'/home/ladu/Code/ENCODE_DREAM_LINK/train_test/train_regions.blacklistfiltered.merged.bed'
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
batch_size = 48


data_input = DreamDataInput.DreamDataLoader(
                                coords_file,
                                None,
                                genome_file,
                                genome_chromosome_sizes,
                                dna_shape_dir,
                                dnase_seq_dir,
                                None,
                                genome_annotation_file,
                                chip_seq_file,
                                dna_window = 600,
                                chrom_window = 3600,
                                chrom_pool_size=8
                                )




seq_len = 600
shape_len = 600
dnase_seq_len = 3600//8
num_classes = 2
num_shapes =4


output_file = 'EGR1_test_batch.h5'

#Pull one batch from eval pool to determine data dimensions
# prior to encoding
(dna_seq_batch,
dna_shape_batch,
dnase_seq_batch,
chip_labels) = data_input.collection.pull_batch_eval(1)

shapes_list = determine_shapes(dna_seq_batch,
                    dna_shape_batch,
                    dnase_seq_batch,
                    chip_labels)


len_1d = shapes_len_1d(shapes_list)
num_examples = data_input.num_examples


with h5py.File(output_file,'w') as hf:
    hf.attrs["dna_seq_len"]= seq_len
    hf.attrs["dna_shape_len"]= shape_len
    hf.attrs["num_dna_shapes"] = num_shapes
    hf.attrs["dnase_seq_len"]= dnase_seq_len
    hf.attrs["num_classes"]= num_classes

    dset = hf.create_dataset('data',
                            (num_examples,len_1d),
                            chunks=(64,len_1d),
                            compression=None,
                            maxshape = (1000000,len_1d))

    #Chunk size of 256 x len_1d is roughly 5 megabytes
    #I found that 64xlen_1d (about 1.25 MB) gives the greatest speed
    #Adding compression can make the read time as much as 5X slower 



    for i in range (num_examples):
    #for i in range(20):
        if i%1000==0:
            print "Wrote",i,"entries"
        #t0 = time.clock()
        (dna_seq_batch,
        dna_shape_batch,
        dnase_seq_batch,
        chip_labels) = data_input.train.pull_batch_train(1)

        shapes_list = determine_shapes(dna_seq_batch,
                        dna_shape_batch,
                        dnase_seq_batch,
                        chip_labels)

        batch_1d = encode_nd_to_1d(dna_seq_batch,
                        dna_shape_batch,
                        dnase_seq_batch,
                        chip_labels)

        dset[i,:]=batch_1d


    #End with
    print "Finished writing data to file", output_file
