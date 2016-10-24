import os
import gzip

#Generate indices for contiguous spans of genomic data defined by files
# ending in "merged.bed"
#And retrieve coordinates by index algorithmically
#Saves on memory and speed when loading genomic coordinates


#If indices exceed the bounds of the total number of elements,
# retrieval should wrap around 

def main():
    input_fname = "ENCODE_DREAM_LINK/train_test/train_regions.blacklistfiltered.merged.bed"
    indexer = IndexCoords(input_fname,200,50)
    
    print indexer.retrieve_by_index(0)
    print indexer.retrieve_by_index(1)
    
    print indexer.retrieve_by_index(5000123)
    print indexer.retrieve_by_index(500000000)
    print indexer.num_elements

    # Raw file detected 51676736 lines
    #Algorithm has 51676481 indices
    #Algorithm is short 255 lines,
    #exactly the same as the number of elements

    
    
class IndexCoords:
    def __init__(self,fname,window_width, step_size):
        self.filename = fname
        self.window_width = window_width
        self.step_size = step_size
        self.intervals = self.load_merged_bed()
        self.num_elements = self.sum_intervals()


        
    #Multiplier will make it possible to 
        
    def load_merged_bed(self):
        intervals = []
        ext= os.path.splitext(self.filename)[1]
        if ext == '.gz' or ext == '.gzip' or ext == '.bgz':
            f =  gzip.open(self.filename,'r')
        elif ext == '.bed':
            f = open(self.filename,'r')
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
                return (interval.contig, start,start+self.window_width,wrap_ind)
        
class MergedCoordsInterval:
    def __init__(self,contig,start,end,window_width,step_size):
        self.contig = contig
        self.start = start
        self.end=end
        self.window_width = window_width
        self.step_size =step_size
        #Note: I never quite figured out why I need to add self.step_size
        nuc_len_less_win = (self.end-self.start-self.window_width)+self.step_size
        self.num_elements = ((nuc_len_less_win)/self.step_size)
        if (nuc_len_less_win % self.step_size != 0):
            print "Error! num_elements must be divisible by step_size"
    
    def retrieve_by_local_index(self,index):
        return self.start+index*self.step_size
    
if __name__ =="__main__":
    main()
