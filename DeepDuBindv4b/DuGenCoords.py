import os,sys
import gzip
#import time
#import dutabix
import h5py
import numpy as np
import random
import glob
import DreamDataInput
import pysam

#Generate indices for contiguous spans of genomic data defined by files
# ending in "merged.bed"
#And retrieve coordinates by index algorithmically
#Saves on memory and speed when loading genomic coordinates


#If indices exceed the bounds of the total number of elements,
# retrieval should wrap around 

def main():
    #merged_bed_fname = "ENCODE_DREAM_LINK/train_test/train_regions.blacklistfiltered.merged.bed"
    #chip_fname = "ENCODE_DREAM_LINK/ChIPseq/labels/ATF7.train.labels.tsv.bgz"
    #indexer = IndexCoords(merged_bed_fname,200,50)
    
    #print indexer.retrieve_by_index(0)
    #print indexer.retrieve_by_index(1)
    #t0 = time.clock()
    #print indexer.retrieve_by_index(5000123)
    #print "Retrieval 1:",(time.clock()-t0)
    #print indexer.retrieve_by_index(900000000)
    #print "Retrieval 2:",(time.clock()-t0)
    #print indexer.num_elements


    #out_dirname = os.path.dirname(chip_fname)
    #index_fname = '.'.join(os.path.basename(chip_fname).split('.')[0:2]+['train_regions.filtered_indices.h5'])

    #print 'Output file is', index_fname
    #culler = IndexCuller(indexer,chip_fname,outfile = index_fname)

    
    
    # Raw file detected 51676736 lines
    #Algorithm has 51676481 indices
    #Algorithm is short 255 lines,
    #exactly the same as the number of elements

    merged_bed_fname = "ENCODE_DREAM_LINK/train_test/train_regions.blacklistfiltered.merged.bed"
    chip_dir = "ENCODE_DREAM_LINK/ChIPseq/labels/"
    indexer = IndexCoords(merged_bed_fname,200,50)
    os.chdir(chip_dir)
    for fname in glob.glob("*.tsv.bgz"):
        output_fname = '.'.join(os.path.basename(fname).split('.')[0:2]+['train_regions.filtered_indices.h5'])
        if os.path.isfile(output_fname):
            print output_fname, " already exists. Skipping index culling"
        else:
            print "Culling indices for", fname
            print 'Output file is', output_fname
            _ = IndexCuller(indexer,fname,outfile = output_fname)
   
       

    
    
class IndexCoords:
    def __init__(self,fname,window_width, step_size):
        self.filename = fname
        self.window_width = window_width
        self.step_size = step_size
        self.intervals = self.load_merged_bed()
        self.num_elements = self.sum_intervals()
        self.iter_ind = 0
       
    #Multiplier will make it possible to 
        
    def load_merged_bed(self):
        intervals = []
        ext= os.path.splitext(self.filename)[1]
        if ext == '.gz' or ext == '.gzip' or ext == '.bgz':
            f =  gzip.open(self.filename,'r')
        elif ext == '.bed':
            f = open(self.filename,'r')
        else:
            print "Merged bed file",self.filename," not found!"
        for line in f:
            contig,start,end = line.strip().split()
            intervals.append( MergedCoordsInterval(contig,
                                                   int(start),
                                                   int(end),
                                                   self.window_width,
                                                   self.step_size) )
        f.close()
        return intervals

    def sum_intervals(self):
        sum_intervals =0
        for interval in self.intervals:
            sum_intervals += interval.num_elements
        return sum_intervals

    def index_from_coords(self, contig, start, end):
        index = 0
        found_index = False
        for interval in self.intervals:
            #print interval.start,interval.end
            if (contig==interval.contig and
                start >= interval.start and
                end <= interval.end):
                dist_to_start = start-interval.start
                if(dist_to_start%interval.step_size!=0):
                    print "Error. start coord should be divisible by step size"
                interval_index = int(dist_to_start//interval.step_size)
                #interval.window_width
                #interval.step_size
                return index+interval_index
            
            #If you passed an interval wihtout
            index += interval.num_elements
       #If you haven't returned a value by this point then there is no index
        #print "Correct index not identified from coordinates" 
        return None    
            
        

    def identify_cell_index(self,index):
        cell_index = 0
        while True:
            block_start = cell_index*self.num_elements 
            block_end = block_start+self.num_elements
            if index >= block_start and index < block_end:
                return_val = cell_index
                break
            else:
                cell_index += 1
            if cell_index >10000:
                print "Error in identify_cell_index"
                print "Must break"
                break
            
        return return_val
    
                     
    def retrieve_by_index(self,index):
        lower =0
        upper=0
        #TODO: add wrap around code
        #If the index number is larger than the number of elements,
        # this may be due to duplicating coords to account for multiple
        # cell types in the training set.
        # In these cases, determine how many wrap arounds have occured
        # (wrap_ind can be used to track cell types).
        
        wrap_ind = 0
        while (index > self.num_elements):  
            index -= self.num_elements
            wrap_ind += 1
        
        for interval in self.intervals:
            lower = upper
            upper += interval.num_elements
            #Determine which interval the current index is on
            if index >= lower and index <= upper:
                start = interval.retrieve_by_local_index(index-lower)
                return (interval.contig, int(start),int(start+self.window_width),int(wrap_ind))

    def __iter__(self):
        return self

    def next(self):
        ret_coords = self.retrieve_by_index(self.iter_ind)
        self.iter_ind += 1
        return ret_coords
            
class MergedCoordsInterval:
    def __init__(self,contig,start,end,window_width,step_size):
        self.contig = contig
        self.start = start
        self.end = end
        self.window_width = window_width
        self.step_size =step_size
        #Note: I never quite figured out why I need to add self.step_size
        nuc_len_less_win = (self.end-self.start-self.window_width +self.step_size)
        self.num_elements = ((nuc_len_less_win)/self.step_size)
        if (nuc_len_less_win % self.step_size != 0):
            print "Error! num_elements must be divisible by step_size"
    
    def retrieve_by_local_index(self,index):
        return self.start+index*self.step_size




class IndexCuller:
    """
    Typically most of the genome is not going to be bound by a particular
    transcription factor. This leads to a heavily skewed data class.
    The labels blacklister should be able to remove the indices of a dataset
    which corresponds to the unbound class. By cutting out
    proportion of the genome that is used for training.
    
    """

    """Visit every entry, by calling the indexer.
    - Record the indices of entries that have B columns and the indicies of entries
       without B columns.
    - Calculate size of B list
    - Randomly remove elements from U list
    - Intercalate B and U list and return result as list of indices.
    - Have option to save this list as a file on a specified directory for
      future use.
    
    """
    
    def __init__(self,indexer,labels_bgz_fname,outfile=None):
        self.indexer = indexer
        self.base_num_elements = self.indexer.num_elements
        self.labels_filename = labels_bgz_fname
        self.U_indices = []
        self.BA_indices = []      
        self.output_file = outfile
        #self.cell_types = DreamDataInput.ChipSeqData.determine_cell_types(self.labels_filename)
        #self.num_cell_types = len(self.cell_types)
        #print "Cell types are," self.cell_types 
        self.num_ba_examples = 0
        self.populate_UB()
        self.cull_U_indices()
        self.write_indices_to_file()
        

    def get_indices(self):
       return  np.concatenate((self.U_indices,self.BA_indices),axis=0)
        
        
    def populate_UB(self):
        #If a row of a chip-seq labels file has only U's, append to U_indices
        # else append to B_indices
        
        #tb = tabix.open(self.labels_filename)
        tb = pysam.TabixFile(self.labels_filename)
        for i in range(self.indexer.num_elements):
        #for i in range(20): 
            diff = 25
            contig,start,_,_ = self.indexer.retrieve_by_index(i)
            #records = tb.query(contig,start-diff,start+diff)
            #records=dutabix.query_popen(self.labels_filename,contig,start-diff,start+diff)
            #sys.stdout.flush()
            records = tb.fetch(contig,start-diff,start+diff)
            if i% 100000 ==0:
                
                print "Visited %d records" % (i)                
            for record in records:
                record = record.strip().split()
                labels = record[3:]
                #print labels
                for j,labl in enumerate(labels):
                    #print j
                    if labl== 'B' or labl =='A':
                        self.BA_indices.append(j*self.base_num_elements+i)
                    else:
                        self.U_indices.append(j*self.base_num_elements+i)

        #print self.U_indices
        #print self.BA_indices

    def cull_U_indices(self):
        #Randomly cull indices for genomic positions which lacks
        # 'A' or 'B' 
        num_ba_examples = len(self.BA_indices)
        num_u_examples = len(self.U_indices)

        if num_u_examples > num_ba_examples:
            #Shuffle in place
            random.shuffle(self.U_indices)
            self.U_indices = self.U_indices[0:num_ba_examples]
        else:
            print "Num BA examples is greater than num U examples"

    def write_indices_to_file(self):
        if self.output_file == None:
            print "Output file not specified"
            print "Not saving indices"
        else:
            print "Saving selected indices to",self.output_file,"!"
            with h5py.File(self.output_file,'w') as hf:
                #
                output_key = "selected_indices" 
                hf.create_dataset(output_key,
                                  data=self.get_indices(),
                                           chunks=True)
        
            
                   
    
            
            
                        
    
if __name__ =="__main__":
    main()
