import numpy as np
import pysam
import os
import DuGenCoords
import glob
import gzip
from DreamDataInput import ChipSeqLabels
import random
import math
import h5py
import time

def main():
    #Purpose: Identify the center of DNAse peaks in training cell types, then check if
    #the center of peaks are bound by tf of interest.
    #Record indices for three categories: bound peaks, unbound peaks, unbound non-peaks
    #Save the indices for the coordinates of these three categories to an h5 file
    #Break up the indices into the following manner: 50% bound peaks, 40% unbound peaks,
    #10% non-peaks

    """
    Order of operation:
    1. Visit every dnase peak relaxed narrowPeak file
    2. Look up center of each peak and check for binding status
    3. If it's bound in more than one cell type, save that index as many times as it shows up
    4. 
    """


    merged_bed_fname = "../ENCODE_DREAM_LINK/train_test/train_regions.blacklistfiltered.merged.bed"
    cell_types = ['GM12878','H1-hESC','HCT116','MCF-7']
    chip_labels_file='../ENCODE_DREAM_LINK/ChIPseq/labels/EGR1.train.labels.tsv.bgz'
    dnase_peaks_dir = '/home/ladu/Code/ENCODE_DREAM_LINK/DNASE/peaks/relaxed/'
    
    indexer = DuGenCoords.IndexCoords(merged_bed_fname,200,50)

    #print indexer.next()
    #print indexer.iter_ind
    #print indexer.next()
    #print indexer.iter_ind
    
    contig,start,end,ind = indexer.retrieve_by_index(12516)
    train_file = 'EGR1.train.indices.h5'
    eval_file = 'EGR1.eval.indices.h5'
    dnase_filt = DnasePeakFilter(indexer,
                                 cell_types,
                                 chip_labels_file,
                                 dnase_peaks_dir,
                                 output_indices_file=train_file,
                                 is_eval_set=False
                                 )
    
    dnase_eval_filt = DnasePeakFilter(indexer,
                                 cell_types,
                                 chip_labels_file,
                                 dnase_peaks_dir,
                                 output_indices_file=eval_file,
                                 is_eval_set=True
                                 )
    

class DnasePeakFilter:
    
    def __init__(self,indexer,cell_types,chip_labels_bgz_file,
                      dnase_peaks_dir,output_indices_file,is_eval_set=False):
        #self.dnase_file_template = ['DNASE','cell_type','relaxed.narrowPeak.sorted.bgz']
        self.dnase_file_template = [ 'DNASE','cell_type','relaxed.narrowPeak']
        self.indexer = indexer

        self.output_indices_file = output_indices_file
        self.output_coords_file = os.path.splitext(output_indices_file)[0]+'.tsv'
        self.chip_labels_file = chip_labels_bgz_file
        f_cell_types = ChipSeqLabels.determine_cell_types(self.chip_labels_file )
        if cell_types != None:
            self.cell_types = cell_types.strip().split(',')
        elif cell_types == None:
            print "(DnasePeakFilter) Cell types not specified."

        if (is_eval_set==True):
            print "(DnasePeakFilter) Evaluation set will use peaks from ChipLabelsFile!"
            #For eval set use all DNASE/ChipLabels data available to select peaks
            self.cell_types = f_cell_types

        print "(DnasePeakFilter) Cell types in",self.chip_labels_file,"are:",f_cell_types
        print "(DnasePeakFilter) Selected cell types are", self.cell_types
        #Locate the proper column 
        self.cell_type_ind_dict = {}
        #print f_cell_types
        for cell_type in self.cell_types:
            #print cell_type
            self.cell_type_ind_dict[cell_type] = f_cell_types.index(cell_type)

                
                
        self.dnase_peaks_dir= dnase_peaks_dir

        self.dnase_peaks_files =  []
        for cell_type in self.cell_types:
            full_dnase_fname = (self.dnase_peaks_dir+os.sep+
                                self.dnase_file_template[0]+'.'+
                                cell_type+'.'+
                                self.dnase_file_template[2])
            print "Examining",full_dnase_fname
            self.dnase_peaks_files.append(full_dnase_fname)

        self.chip_tabix = pysam.TabixFile(self.chip_labels_file)
        self.bound_peak_indices = []
        self.unbound_peak_indices = []
        self.unbound_nonpeak_indices =[]

        
        if is_eval_set ==False:
            self.locate_bound_peaks()
            self.shuffle_normalize_indices()
            self.write_indices_to_hdf5(self.get_shuffled_indices())
            #self.write_coords_to_tsv(self.get_shuffled_indices())
        elif is_eval_set == True:
            self.locate_dnase_peaks_only(include_random_coords=True)
            print "Evaluation set selected"
            print "Selecting genomic locations surrounding dnase peaks"
            #For evaluating datasets (such as ladder or final submission data),
            #we need to evaluate the genomic region surrounding each peak

            #The eval window is defined in indices. If the step size is 50 bp, then
            # a window [-10,12] would represent downstream 500 and upstream 600
            #
            #Note in order to get a submission in on time, I had to make this -2,+2
            self.eval_window = [-2,2]
            eval_indices = self.get_indices_window(self.eval_window,collapse_cell_types=True)
            self.write_indices_to_hdf5(eval_indices)
            #self.write_coords_to_tsv(eval_indices)

            
    def locate_bound_peaks(self):
        #Compare ChIP labels file and dnase relaxed peaks file and identified bound/unbound state
        # Using an IndexCoords object, identify the indices
        
        for cell_type in self.cell_types:
            f = (self.dnase_peaks_dir+os.sep+
                            self.dnase_file_template[0]+'.'+
                            cell_type+'.'+
                            self.dnase_file_template[2])


            ext =os.path.splitext(f)[1]
            if ext == '.bgz' or ext == '.gz' or ext =='.gzip':
                fhandle = gzip.open(f,'r')
            else:
                fhandle = open(f,'r')


            #Prev end is used to pick a random spot between features
            #This is used for populating the unbound_nonpeaks list
            prev_end = 0
            #line_count=0
            while True:
                
                #line_count += 1
                line = fhandle.readline()
                if line == "":
                    print "EOF reached"
                    break

                line = line.strip().split()
                point_peak = int(int(line[1])+int(line[9]))


                #Round start and end to nearest 50 surrounding point peak
                contig = line[0]
                start = int(50*(point_peak//50))
                end = start+50


                cell_ind = self.cell_type_ind_dict[cell_type]
                chip_label_col = 3+cell_ind
                
                train_index = self.indexer.index_from_coords(contig,start,end)
                #print (contig,start,end,train_index)
                try:
                    chip_records = self.chip_tabix.fetch(contig,start-5,start+5)
                    rec = chip_records.next().strip().split()
                    #print "Record for",contig,start-5,start+5,"located"
                except (StopIteration,ValueError):
                    rec = None

                
                if rec!= None and (rec[chip_label_col] == 'B' or rec[chip_label_col]=='A'):
                    #print 'Bound found!'
                    #Note need to handle coords not present
                    if train_index == None:
                        #print "Index for coords",contig,start,end, "not found."
                        pass
                    else:
                        
                        self.bound_peak_indices.append(train_index+
                                                       self.indexer.num_elements*cell_ind)
                elif rec !=None and rec[chip_label_col] == 'U':
                    if train_index == None:
                        #print "Index for coords",contig,start,end, "not found."
                        pass
                    else:
                        #print 'U'
                        self.unbound_peak_indices.append(train_index+
                                                         self.indexer.num_elements*cell_ind)
                #Find a random spot between prev_end and start and pick that location
                # for populating unbound_nonpeaks
                if start-prev_end > 400:
                    random_location = random.randrange(prev_end,start)
                    rand_start = int(50*(random_location//50))
                    rand_end = rand_start+50
                    
                    try:
                        rand_records = self.chip_tabix.fetch(contig,start,end) 
                        rand_rec = rand_records.next().strip().split()
                    except (StopIteration,ValueError):
                        rand_rec = None


                    if rand_rec != None and rand_rec[chip_label_col] == 'U':
                        crap_index = self.indexer.index_from_coords(contig,rand_start,rand_end)
                        if crap_index == None:
                            #print "Index for crappy coords",contig,rand_start,rand_end, "not found."
                            pass
                        else:    
                            self.unbound_nonpeak_indices.append(crap_index+
                                                             self.indexer.num_elements*cell_ind)
                prev_end = end



    def locate_dnase_peaks_only(self,include_random_coords=True):
        for cell_type in self.cell_types:
            f = (self.dnase_peaks_dir+os.sep+
                            self.dnase_file_template[0]+'.'+
                            cell_type+'.'+
                            self.dnase_file_template[2])
            ext =os.path.splitext(f)[1]
            if ext == '.bgz' or ext == '.gz' or ext =='.gzip':
                fhandle = gzip.open(f,'r')
            else:
                fhandle = open(f,'r')
            #Prev end is used to pick a random spot between features
            #This is used for populating the unbound_nonpeaks list
            prev_end = 0
            #line_count=0
            while True:
                
                #line_count += 1
                line = fhandle.readline()
                if line == "":
                    print "EOF reached"
                    break

                line = line.strip().split()
                point_peak = int(int(line[1])+int(line[9]))


                #Round start and end to nearest 50 surrounding point peak
                contig = line[0]
                start = int(50*(point_peak//50))
                end = start+50

                #There is no labels file. Just save all DNAse peaks
                cell_ind = self.cell_type_ind_dict[cell_type]
                
                eval_index = self.indexer.index_from_coords(contig,start,end)
                if eval_index == None:
                    #Your in a blacklisted region of the genome, do not bother here
                    continue
                

                self.unbound_peak_indices.append(eval_index+
                                                         self.indexer.num_elements*cell_ind)

                #Find a random spot between prev_end and start and pick that location
                # for populating unbound_nonpeaks

                if start-prev_end > 400:
                    random_location = random.randrange(prev_end,start)
                    rand_start = int(50*(random_location//50))
                    rand_end = rand_start+50

                    rand_index = self.indexer.index_from_coords(contig,rand_start,rand_end)
                    if rand_index == None:
                    #Your in a blacklisted region of the genome, do not bother here
                        continue
                    self.unbound_peak_indices.append(rand_index+
                                                         self.indexer.num_elements*cell_ind)
                prev_end = end

                
    def shuffle_normalize_indices(self):
        #50% bound peaks, %40 unbound peaks, 10% unbound random locations
        #Always keep all bound peaks
        frac_unbound_peaks =.4
        frac_unbound_nonpeaks = .1
        self.bound_peak_indices  = np.random.permutation(self.bound_peak_indices)
        self.unbound_peak_indices = np.random.permutation(self.unbound_peak_indices)
        self.unbound_nonpeak_indices = np.random.permutation(self.unbound_nonpeak_indices)           

        num_bound_peaks = len(self.bound_peak_indices)
        num_unbound_peaks = len(self.unbound_peak_indices)
        num_unbound_nonpeaks = len(self.unbound_nonpeak_indices)
        print "Total number of bound peaks identified:",num_bound_peaks
        print "Total number of unbound peaks identified:",num_unbound_peaks
        print "Total number of unbound non peaks identified:",num_unbound_nonpeaks
        
        if (num_bound_peaks<=num_unbound_peaks):
            num_unbound_peaks = int(math.floor(num_bound_peaks*frac_unbound_peaks*2))
            num_unbound_nonpeaks = int(math.floor(num_bound_peaks*frac_unbound_nonpeaks*2))
            print "Number of selected bound peaks:",num_bound_peaks
            print "Number of selected unbound peaks:",num_unbound_peaks
            print "Number of selected unbound non-peaks:",num_unbound_nonpeaks
            self.unbound_peak_indices = self.unbound_peak_indices[0:num_unbound_peaks]
            self.unbound_nonpeak_indices = self.unbound_nonpeak_indices[0:num_unbound_nonpeaks]
        elif (num_bound_peaks>num_unbound_peaks):
            print "Number of bound peaks is greater than number of unbound peaks"
            num_unbound_non_peaks = num_bound_peaks-num_unbound_peaks
            self.unbound_nonpeak_indices=self.unbound_nonpeak_indices[0:num_unbound_nonpeaks]
            print "Number of selected bound peaks:",num_bound_peaks
            print "Number of selected unbound peaks:",num_unbound_peaks
            print "Number of selected unbound non-peaks:",num_unbound_nonpeaks


    def get_unshuffled_indices(self):
        #Actually each of the three categories being concatenate are
        #  already shuffled internally
        all_indices = np.concatenate((self.unbound_nonpeak_indices,
                                      self.bound_peak_indices,
                                      self.unbound_peak_indices),axis=0)
        return all_indices    

        
    def get_shuffled_indices(self):
        all_indices = np.concatenate((self.unbound_nonpeak_indices,
                                      self.bound_peak_indices,
                                      self.unbound_peak_indices),axis=0)
        all_indices = np.random.permutation(all_indices)
        return all_indices

        
        
    
    def get_indices_window(self,eval_window,collapse_cell_types = False):
        t0 = time.clock()
        #np.unique also sorts
        indices = np.unique(self.get_unshuffled_indices())
        #Creating evaluation windows
        #Preallocate array
        surround_size = -eval_window[0]+eval_window[1]
        new_indices = np.zeros(surround_size*len(indices))

        
        for ind,i in enumerate(indices):
            if collapse_cell_types == True:
                if (i > self.indexer.num_elements):
                    #Remap celltype
                    cell_index = self.indexer.identify_cell_index(i)
                    #Subtract overflor from index 
                    i = i-(cell_index*self.indexer.num_elements)
            
            surrounding_indices= np.arange(eval_window[0]+i,
                                           eval_window[1]+i)
            new_indices[(ind*surround_size):((ind+1)*surround_size)] =surrounding_indices
        new_indices = np.unique(new_indices)
        print "get_indices_window took", time.clock()-t0, "s to run"
        print "(indices window for eval) Num elements:",len(new_indices)
        return new_indices.clip(0)

                
    def write_indices_to_hdf5(self,indices):
        if (self.output_indices_file == None) or (self.output_indices_file==''):
            print "Output file not specified"
            print "Not saving indices"
        else:
            print "Saving selected indices to",self.output_indices_file,"!"
            with h5py.File(self.output_indices_file,'w') as hf:
                output_key = "selected_indices" 
                hf.create_dataset(output_key,
                                data=indices,
                                       chunks=True)
        
    def write_coords_to_tsv(self,indices):
        print "Saving coordinates to",self.output_coords_file,"!"
        indices = np.unique(indices)
        with open (self.output_coords_file,'w') as f:
            for i in indices:
                tup = self.indexer.retrieve_by_index(i)
                tup = '\t'.join([tup[0],str(tup[1]),str(tup[2]),str(tup[3]),'\n'])
                f.write(tup)
        
if __name__ == "__main__":
    main()

   
