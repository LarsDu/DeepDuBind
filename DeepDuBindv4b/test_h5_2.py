import numpy as np
import h5py
import os
import time
import DreamDataSaverHdf5
import DreamDataInput
import gc

#test_fname = '../ENCODE_DREAM_DATASETS/EGR1_test_batch.h5'
test_fname = './REST_train_infC/REST.train_cache.h5'

data_name='data'





with h5py.File(test_fname,'r') as hf:

    print hf.attrs.get("seq_len")
    print hf.attrs["dna_shape_len"]
    print hf.attrs["pooled_chrom_window"]
    print hf.attrs["num_classes"]

    dset = hf[data_name]
    shape_list = [(1,1,hf.attrs["seq_len"],4),
                  (1,1,hf.attrs["dna_shape_len"],4),
                  (1,hf.attrs["pooled_chrom_window"]),
                  (1,hf.attrs["num_classes"])]
    

    print dset.chunks
    print dset.shape
    
    t0 = time.clock()
    print  dset[0]
    print "Elapsed:",time.clock()-t0

    t0 = time.clock()
    print  dset[18]
    print "Elapsed:",time.clock()-t0


    t0 = time.clock()
    print  dset[102401]
    print "Elapsed:",time.clock()-t0


    
    t0 = time.clock()
    print  dset[200001]
    print "Elapsed:",time.clock()-t0

    
    t0 = time.clock()
    vals =DreamDataSaverHdf5.decode_1d_to_nd(dset[200200],shape_list)
    print "pulled stuff"
    print "Elapsed:",time.clock()-t0

    t0 = time.clock()
    vals =DreamDataSaverHdf5.decode_1d_to_nd(dset[277080],shape_list)
    print "pulled stuff"
    print "Elapsed:",time.clock()-t0

    print vals[0].shape,vals[1].shape,vals[2].shape,vals[3].shape
    print np.reshape(vals[2],(1,1,450,1)).shape

    print vals[0]
    print vals[1]
    print vals[2]
    print vals[3]

    t0 = time.clock()
    randnums = np.random.random_integers(6400,360000,64)
    for i in randnums:
        vals = DreamDataSaverHdf5.decode_1d_to_nd(dset[i],shape_list)
        #print vals[3]
                                                      
    t1= time.clock()

    print "64 Records time:",t1-t0

    
training_set = DreamDataInput.DataCollectionHdf5(test_fname)
training_set.open()
start_time = time.time()

training_set.pull_batch_train(64)
#feed_dict = {"dna_seq_placeholder":training_set.dna_seq_batch,
#            "dna_shape_placeholder":training_set.dna_shape_batch,
#            "dnase_seq_placeholder":training_set.dnase_seq_batch,
#             "chip_labels_placeholder":training_set.chip_labels_batch} 

#training_set.pull_batch_train(64)
duration = time.time()-start_time
print "Time taken for batch retrieval", duration

#t0 = time.time()
#training_set.pull_batch_train(64)
#print "Time for just pulling a batch without save",time.time()-t0
training_set.close()



