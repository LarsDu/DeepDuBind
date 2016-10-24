from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import motifs
import gzip

import numpy as np
import random
#import multiprocessing 
import itertools

import os
import sys
"""
A class for holding static methods that can serve as useful
Bioinformatics tools. There are methods here for random sequence
generation and tools for converting sequence into one-hot binary
vectors for machine learning applications.     


"""

def seq_permutation(nuc_len):
    nucs = 'ATGC'*nuc_len
    all_perm = itertools.permutations(nucs,nuc_len)
    for i in all_perm:
        print ''.join(i)



def get_str_from_coords(SeqObj,start,end):
    #print SeqObj
    return str(SeqObj[start:end])


def seq_to_onehot_topdown(seqObj):
    """
    Converts a Seq object into a set of four boolean one-hot vectors
    vectorDict order is A T G C. Remember the proper command
    for transposing a numpy array is myarray.transpose()
    """
    vectorDict = {'A': [1, 0, 0, 0],
                  'a': [1, 0, 0, 0],
                  'T': [0, 1, 0, 0],
                  't': [0, 1, 0, 0],
                  'G': [0, 0, 1, 0],
                  'g': [0, 0, 1, 0],
                  'C': [0, 0, 0, 1],
                  'c': [1, 0, 0, 1],
                  'N': [0.25, 0.25, 0.25, 0.25],
                  'n': [0.25, 0.25, 0.25, 0.25]
                  }

    seq_str = str(seqObj)
    onehot = np.zeros((len(seq_str),4),dtype='float32')
    for i,letter in enumerate(seq_str):
        onehot[i,:] = np.array(vectorDict[letter])
    return onehot


def seq_to_onehot(seqObj):
    #Converts a sequence to a 4 one-hot vectors in the order ATGC
    seq = seq_to_onehot_topdown(seqObj)
    return seq.transpose()

###Methods for generating random nucleotide sequences, and random nucleotide
### sequences seeded with user specified motifs. These will be used to generate
### simulated training sets for validating machine learning efficacy


def seq_to_flat_onehot(seqObj):
    #Does the same thing as seq_to_onehot, but flattens the output
    seq = seq_to_onehot(seqObj)
    return seq.ravel()


@staticmethod
def fasta_to_flat_onehot(fname):
# Takes in fasta file and converts to numpy array
# > First output is a large numpy array with each
#   row representing a flattened one-hot vector.
#   Nucleotides are in the order A,T,G,C
# Usage:
#  fasta_to_flat_onehot('myfasta.fa')


    seq_parser = SeqIO.parse(fname,"fasta")

    #Convert generator fasta objects to large list of
    #1xn (n=4*SEQLEN) onehot vectors
    seq = seqio_to_flat_onehot(seq_parser)

    #Close iterator handles
    seq_parser.close()

    return seq



def seqio_to_flat_onehot(bio_seqrecord_gener):
# Takes a BioPython SeqIO generator, converts every element into a
# 4 rows of one-hot vectors, with each row going up to down being A, T,
# G, then C. Then flatten/reshape, into 1x(n*4) numpy array vector.
#The second return value is the length of each nucleotide string
#The flat_onehot return value is a numpy array with each row being a different
#flattened one-hot vector representing a single fasta record.

    #Examine first element to determine sequence length
    rec0seq = bio_seqrecord_gener.next().seq
    nuc_len = len(rec0seq)
    vec_len = nuc_len*4

    #Initialize width of output 
    flat_onehot = np.zeros((1,vec_len), dtype=np.bool)
    #Set first row
    flat_onehot[0,:] =seq_to_flat_onehot(rec0seq)


    #Set the remaining records in the generator
    for rec in bio_seqrecord_gener:
        #Convert to one-hot
        cur_rec = seq_to_flat_onehot(rec.seq)
        #Flatten and append as row to 
        flat_onehot = np.vstack((flat_onehot,cur_rec))

    return flat_onehot





def rand_dna_nuc(nuc_len, gc_fraction=0.5):
    # Generates a random nucleotide sequence of fixed length with 
    # specified nucleotide percentage
    # requires numpy

    s_perc= (gc_fraction*.5)
    w_perc = (.5*(1-gc_fraction))
    nucs = ['A','T','G','C']
    probs = [s_perc,s_perc,w_perc,w_perc]

    return ''.join([np.random.choice(nucs,p=probs)
                     for _ in xrange(nuc_len)])




def shuffle_batch_pull(batch_size,*args):
# Takes in a user specified batch_size number,
# and one or more numpy arrays that have the same number of rows
# Randomly pulls batch_size number of rows from input numpy_arrays
# args* here can be a variable number of numpy arrays
# This function will return random rows (with the same index)
# from each input numpy array.
# If you want to pull data from a different axis, you can specify that axis

#This function is used for pulling random batches for stochastic gradient
#descent

   #Check if there are the same number of rows in each input numpy array
    AXIS = 0
    num_rows = args[0].shape[AXIS]
    random_rows = np.random.choice(num_rows,batch_size, replace=False)

    return_tuple = []
    for i,np_arr in enumerate(args):
        if np_arr.shape[0]!=num_rows:
            print 'Numpy array row mismatch error. Exiting!'
            return None
        else:
            return_tuple.append( np.take(np_arr, random_rows,axis=AXIS))

    return tuple(return_tuple)






def rand_motif_instance(motif):
    #Takes a motif object, and pulls an instance randomly from it
    #Note: A motif object can be created by 
    #mymotif = motifs.create(instances)
    #where instances = [Seq("ATA"),Seq("AGA")]
    return random.choice(motif.instances)


def rand_seed_motif(nuc_seq, motif, dist_from_tss=[-50,50],
                    orientation= 0,wobble=0):
        #Places a specified motif within a specified nucleotide sequence
        #Optional parameters are:
        #     >dist_from_tss - specifies how far away from center of nuc_seq 
        #     the motif should be placed. Can be a list of positive and negative
        #     values
        #     >orientation = -1 to place reverse orientation, 1 for forward, 0 for 
        #     random orientation              
        #     >wobble = number of nucleotides the motif can be out of register
        #     with respect to dist_from_tss
        #Note: The 5' end of the motif is always inserted at dist_from_tss
        #Regardless of orientation

        #motif needs to be a Biopython motif object

        nuc_len = len(nuc_seq)
        motif_len = len(motif) 
        tss_center = np.floor(nuc_len*.5)


        #Set motif orientation
        if orientation ==0:
                orientation = int(random.choice([-1,1])  )
        if orientation == -1:
                motif = motif.reverse_complement()       
        elif orientation == 1:
                pass
        elif (orientation >1 or orientation < -1):
                print ("Orientation out of bounds. Needs to be -1,0,or 1")
                return None

        insertion_points = [tss_center + each_dist for each_dist in dist_from_tss]
        #Apply wobble to insertion points
        insertion_points = [int(each_point) +random.randint(0,wobble)for each_point  in insertion_points] 
        #Check bounds, exit if dist from tss too far              
        for each_point in insertion_points:
            if each_point < 0 or each_point+motif_len>nuc_len:
                print("Motif insertion point out of bounds. Exiting...")
                return None
            #Insert motif sequence at appropriate position, replacing original sequence        nuc_seq[tss_center-dist_from_tss]
            #Pick a random motif from all instances



            nuc_seq =( nuc_seq[:int(each_point)]+
                       str(motif)+nuc_seq[int(each_point+motif_len):]  )
        return nuc_seq
        #Important python note: my_list[random.randint(1,len(my_list))-1] for
        # randomly selecting an item in a list.
        #Or just use random.choice(my_list)


def motif_to_convfilter(motifObj,filter_width):
    #Return a python list of numpy convolution filters
    #This method will go into all the instances in the bioPython
    # motif, and make a corresponding [1,motif_len,5,1] convolutional
    # filter where dim_2 = 5 represents 5 different nucleotide
    # letters (equivalent to colors in image processing)
    # I am using this method strictly for testing tensorflow's
    # convolution functions on nucleotide sequences

    filter_list = []
    for instance in motifObj.instances:
        filter_list.append(seq_to_convfilter(instance,filter_width))
    return filter_list   


def seq_to_convfilter(seqObj,filter_width=9):
   
    onehot = seq_to_onehot_topdown(seqObj)
    #We want to convert this one_hot of shape [filter_width,4]
    # to a [1,filter_width,4,1] filter (remember the last dim,
    #dim_3 is used by tensorflow in the case you have multiple filters).
    #Use numpy indexing to put this one_hot representation in the correct
    #dims of the 4d matrix
    seq_len = len(seqObj)
    if seq_len>filter_width:
        print "Motif is longer than specified filter width"
        print "quitting"
        return
    elif seq_len == filter_width:
        pass
    elif seq_len<filter_width:
        #If the filter_width>seq_len
        #pad the onehot array with zeros on the top and bottom
        #So that the number of rows == filter_width
        diff = filter_width-seq_len

        #Calculate row index for centering sequence on filter 
        centering_ind = (np.floor(filter_width/2)-np.floor(seq_len/2))
        temp = np.zeros(shape=(filter_width,4))
        #This confusing line of code essentially drops the one_hot motif
        #representation into the middle of the filter
        temp[centering_ind:(onehot.shape[0]+centering_ind),:] = onehot
        onehot=temp
    #Initialize empty ndarray 
    convfilter = np.empty(shape=(1,filter_width,4,1))
    convfilter[0,:,:,0]=onehot
    return convfilter #shape [1,filter_width,4,1]





def seq_to_4d_onehot(seqObj):
    # Convert Seq object (BioPython format) to a numpy ndarray
    # with the shape [batch_size=1,height=1,width=seq_len,num_channels=4]
    onehot = seq_to_onehot_topdown(seqObj) #shape is [seq_len,4]
    return onehot[np.newaxis,np.newaxis,:,:]


def fasta_to_4d_onehot(fname):
    '''
    Converts a fasta file into an ndarray with the shape
    [num_fasta_entries,height=1,seq_len,num_channels =4]

    dim3 (num_channels) is for each nucleotide letter, with the
    ordering A,T,G,C        
    '''
    seq_parser = SeqIO.parse(fname,"fasta")
    seq_dict =SeqIO.index(fname,"fasta")
    num_records = len(seq_dict)
    #Convert generator fasta objects to large list of
    #1xn (n=4*SEQLEN) onehot vectors
    seq = seqio_to_4d_onehot(seq_parser,num_records)

    #Close iterator handles
    seq_parser.close()
    seq_dict.close()
    return seq


def seqio_to_4d_onehot(bio_seqrecord_gener,num_records):
    '''
    Converts a fasta file into an ndarray with the shape
    [num_fasta_entries,height=1,seq_len,num_channels =4]

    dim3 (num_channels) is for each nucleotide letter, with the
    ordering A,T,G,C,N

    num_records parameter is necessary for numpy array preallocation
    '''

    #Use BioPython indexing to determine the number of records in
    #The generator object
    #Examine first element to determine sequence length and
    #use that information to initialize empty numpy ndarray
    first_rec = bio_seqrecord_gener.next().seq
    seq_len = len(first_rec)
    #Preallocate array (this makes things memory efficient)

    test = np.zeros((num_records,1,seq_len,4),dtype='float32')
    onehot_4d = np.zeros((num_records,1,seq_len,4),dtype='float32')
    onehot_4d[0,:,:,:] = seq_to_4d_onehot(first_rec)

    #Set the remaining records in the generator
    for i,rec in enumerate(bio_seqrecord_gener):
        
        #Look at first entry to determine nucleotide sequence length
        #Convert each sequence to one-hot
        onehot_4d[i+1,:,:,:]=seq_to_4d_onehot(rec.seq) 
    return onehot_4d



def onehot_4d_to_nuc(onehot_4d,output_file=None,include_fasta_header=False):
    #Converts ndarray with shape
    #[num_fasta_entries,height=1,seq_len,num_channels =4]
    #into a fasta sequence.

    stdout_old = sys.stdout
    if (output_file !=None):
        sys.stdout = open(output_file,'w')
        
    onehot_4d = np.squeeze(onehot_4d,axis=1)
    #[num_fasta_entries,seq_len,num_channels =4]
    num_entries = onehot_4d.shape[0]
    all_nucs = []
    for i in range(num_entries):
        nuc_string = np.repeat('',onehot_4d.shape[1])    
        A_mask = onehot_4d[i,:,0] == 1
        T_mask = onehot_4d[i,:,1] == 1
        G_mask = onehot_4d[i,:,2] == 1
        C_mask = onehot_4d[i,:,3] == 1
        N_mask = ((onehot_4d[i,:,0] == .25)+
                  (onehot_4d[i,:,1] == .25)+
                  (onehot_4d[i,:,2] == .25)+
                  (onehot_4d[i,:,3] == .25))
        nuc_string[A_mask]='A'
        nuc_string[T_mask]='T'
        nuc_string[G_mask]='G'
        nuc_string[C_mask]='C'
        nuc_string[N_mask]='N'
        if (include_fasta_header == True):
            print '>seq',i
        nucleotides = nuc_string.tostring()
        print nucleotides
        all_nucs.append(nucleotides)
    sys.stdout = stdout_old
    return all_nucs



def extract_n_classes_fasta(fname_list):
    '''
    Convert list of fasta files to one-hot vectors,
    and creates one-hot labels for each class.
    Labels are ordered based on file inputs.
    Each fasta file gets a unique label.
    
    '''
    if type(fname_list) is not list:
        fname_list = [fname_list]    
    rec_onehot = []
    labels_list = []
    num_classes= len(fname_list)

    for i,fname in enumerate(fname_list):
        # fasta_to_4d_onehot Converts each fasta file into an
        # ndarray with the shape:
        # [num_fasta_entries,height=1,seq_len,num_channels =4]
        rec_onehot.append(fasta_to_4d_onehot(fname))
        label_block = np.zeros((rec_onehot[i].shape[0],num_classes),dtype=np.bool_)
        #Convert column corresponding to current index to all ones
        label_block[:,i] = np.ones(label_block.shape[0],dtype=np.bool_)
        labels_list.append(label_block)

    #Vertical stack data
    #Reshape python one-hot list into numpy matrix by vertical stacking of entries.
    data = np.concatenate(rec_onehot,axis=0) 
    labels = np.concatenate(labels_list,axis=0)
    return(data,labels,num_classes)


def divide_fasta(fname_list, fractions = [.6,.2,.2], suffix_list=['train','test','validation']):
    '''
    Divide a fasta file into multiple fasta files with 
    '''
    if len(fractions)!= len(suffix_list):
        print('Size mismatch between fraction division list and suffix list!')
    num_outputs = len(fractions)    

    if sum(fractions)> 1:
        print 'Error! Values in \'fractions\' must add up to <= 1'
        return None
    
    #Validate list
    if type(fname_list) is not list:
        print "Converting to list"
        fname_list = [fname_list]

    for i,file in enumerate(fname_list):
        print 'Opening file', file
        handle = open(file, "rU")
        records = list(SeqIO.parse(handle,"fasta"))
        handle.close()
        random.shuffle(records)
        num_records = len(records)
        print ('File ',file,' contains ', len(records), ' records.')

        prev_record_index = 0
        for k in range(num_outputs):
        #Divide the original data set into the values specified by
        #the fractions list    
            num_cur_records= int(math.floor(num_records*fraction[k])) 
            lower_range = prev_record_index
            upper_range = prev_record_index+num_cur_records
            prev_record_index = upper_range
            _write_fasta_file(file,file+'_'+suffix_list[i],
                        records[lower_range:upper_range])

def retrieve_files_on_dir (dir,ext):
    '''
    Retrive all files on dir with given extension
    '''
    dir_files = os.listdir(dir)
    return [(dir+os.sep+f) for f in dir_files if f.endswith(ext)]
                      
    

            
def _write_fasta_file(input_fname,name_extension,seqio_records_list):
    fname, file_extension = os.path.splitext(input_fname)
    output_fname = fname+'_'+name_extension+file_extension
    output_handle = open(output_fname,"w")
    SeqIO.write(seqio_records_list,output_handle,"fasta")
    output_handle.close()


def count_lines_file(filename):
    '''Counts the number of lines in a file'''
    #From http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    file_ext = os.path.splitext(filename)[1]
    if file_ext == '.gz' or file_ext=='.gzip':
        f = gzip.open(filename,'r')
    else:
        f = open(filename,'r')                  
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization
    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    f.close()
    return lines


def extract_seq_by_coords(genome_fasta, chr_key,coord_tuple):
    """"
    Extract sequence from fasta file. User must provide header key
    (typically the chromosome index). This method opens and loads the file on
    each call, so use it sparingly!
    """
    genome_dict = SeqIO.index(genome_fasta,"fasta")
    return genome_dict[str(chr_key )].seq[coord_tuple[0]:coord_tuple[1]]
    genome_dict.close()


def split_fasta(fasta_file):
    '''Splits a single fasta file into multiple fasta files
       with each output file representing one fasta entry.
       Useful for splitting genomic data into constituent chromosomes'''
    
    for entry in SeqIO.parse(fasta_file,'fasta'):
        output_fname = entry.id+'.fa'
        with open(output_fname,'w') as out_file:
            out_file.write('>'+entry.id+'\n')
            out_file.write(str(entry.seq))



def contig_sizes_to_dict(sizes_fname):
    contig_size_dict = {}
    with open(sizes_fname,'r') as f:
        for line in f:
            contig,contig_len = line.strip().split()
            contig_size_dict[contig]=int(contig_len)

    return contig_size_dict

            
