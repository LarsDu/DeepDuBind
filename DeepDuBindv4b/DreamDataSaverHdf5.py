import DreamDataInput #Pull data batches
import DuJsonParse
import numpy as np
import os
import sys
import time
import h5py


def main():
    pass
        
def create_hdf5(json_params_file,peaks_index_file, output_file):
    print "Creating", output_file
    print "Creating DreamDataLoader object to pull training examples!"
    print "This file may take over an hour to generate!"
    data_input = DreamDataInput.DreamDataLoaderJson(json_params_file,peaks_index_file)
    json_params = data_input.params
    print "Saving to output file", output_file

    #Pull one batch from eval pool to determine data dimensions
    # prior to encoding
    (dna_seq_batch,
    dna_shape_batch,
    dnase_seq_batch,
    chip_labels,
    _) = data_input.collection.pull_batch_eval(1)

    shapes_list = determine_shapes(dna_seq_batch,
                        dna_shape_batch,
                        dnase_seq_batch,
                        chip_labels)
        
    
    len_1d = shapes_len_1d(shapes_list)
   

    
    with h5py.File(output_file,'w') as hf:
        hf.attrs["seq_len"]= json_params.seq_len
        hf.attrs["dna_shape_len"]= json_params.seq_len
        hf.attrs["num_dna_shapes"] = 4
        hf.attrs["pooled_chrom_window"]= json_params.pooled_chrom_window
        hf.attrs["num_classes"]= data_input.num_classes
       
        dset = hf.create_dataset('data',
                                (data_input.num_examples,len_1d),
                                chunks=(64,len_1d),
                                compression=None,
                                maxshape = (1000000,len_1d))

        #Chunk size of 256 x len_1d is roughly 5 megabytes
        #I found that 64xlen_1d (about 1.25 MB) gives the greatest speed
        #Adding compression can make the read time as much as 5X slower 

        

        for i in range (data_input.num_examples):
        #for i in range(20):
            if i%1000==0:
                print "Wrote",i,"entries"
            #t0 = time.clock()
            (dna_seq_batch,
            dna_shape_batch,
            dnase_seq_batch,
            chip_labels, _) = data_input.collection.pull_batch_train(1)

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
    
        
            
def determine_shapes(*argv):
    #Determine the shapes of numpy arrays
    shape_list = []
    for arg in argv:
        shape_list.append(arg.shape)
    return shape_list

def shapes_len_1d(shapes):
    total_len=0
    for shape in shapes:
        total_len += np.prod(shape)
    return total_len

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
    #t0 = time.clock()
    for shape in shape_list:
        elem_len= np.prod(shape)
        return_list.append(np.reshape(input_1d[start_ind:elem_len+start_ind],shape))
        start_ind += elem_len

    #print "Time", time.clock()-t0
    return return_list
        



if __name__=="__main__":
    main()
    
