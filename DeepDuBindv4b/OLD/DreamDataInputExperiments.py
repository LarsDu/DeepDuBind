import numpy as np
import os
import DuBioTools as dbt
#from pybedtools import BedTool
import gzip
import pyBigWig
import pysam
from Bio import SeqIO

import h5py #for dnashape data
import pyBigWig # for bigwig data
#import tabix 
#import dutabix

import time
import DuGenCoords
import DuBioConverter
import DreamDataSaverHdf5
#from DnasePeakFilter import DnasePeakFilter


import DuJsonParse
import gc
import cProfile

class DataCollectionHdf5:
    #Reads preprocessed files generated by DreamDataSaverHdf5
    def __init__(self, hdf5_fname):
        self.filename = hdf5_fname
        self.data_key = 'data'
        
        with h5py.File(self.filename,'r') as hf:
            self.seq_len = hf.attrs["seq_len"]
            self.dna_shape_len = hf.attrs["dna_shape_len"]
            #self.num_dna_shapes = hf.attrs["num_dna_shapes"]
            self.num_dna_shapes = 4
            self.pooled_chrom_window =  hf.attrs["pooled_chrom_window"]
            self.num_classes = hf.attrs["num_classes"]


        self.shape_list= [(1,1,self.seq_len,4),
                          (1,1,self.dna_shape_len,self.num_dna_shapes),
                          (1,self.pooled_chrom_window),
                          (1,self.num_classes)]

        self.fhandle = h5py.File(self.filename,'r')
        self.data = self.fhandle[self.data_key]
        self.num_examples = self.data.shape[0]
        self.perm_indices = np.random.permutation(range(self.num_examples))
        self.num_examples = len(self.perm_indices)
        self.epoch_tracker = EpochTracker(self.num_examples)
        self.eval_epoch_tracker = EpochTracker(self.num_examples)
        

        #To set these values, the file must be opened
        self.fhandle = None
        self.data = None

        self.dna_seq_batch = np.zeros( (64,1,self.seq_len,4) )
        self.dna_shape_batch = np.zeros((64,
                                    1,
                                    self.dna_shape_len,                                    
                                    self.num_dna_shapes) )
        self.dnase_seq_batch = np.zeros( (64,self.pooled_chrom_window) )
        self.chip_labels_batch = np.zeros( (64,self.num_classes) )
        self.batch_indices = np.zeros(64)

    def set_perm_indices(self,new_indices):
        self.perm_indices = new_indices
        self.num_examples = len(self.perm_indices)
        self.epoch_tracker = EpochTracker(self.num_examples)
        self.eval_epoch_tracker = EpochTracker(self.num_examples)
        

        
    def open(self):
        print "Opening HDF5 data file:",self.filename
        self.fhandle = h5py.File(self.filename,'r')
        self.data = self.fhandle[self.data_key]
        
    def close(self):
        "Always remember to call this function when finished!!!"
        self.fhandle.close()


    def _pull_batch(self,batch_size,epoch_tracker):
        t0 = time.clock()
        #Preallocate
        
        #dna_seq_batch = np.zeros( (batch_size,1,self.seq_len,4) )
        #dna_shape_batch = np.zeros((batch_size,
        #                            1,
        #                            self.dna_shape_len,                                    
        #                            self.num_dna_shapes) )
        #dnase_seq_batch = np.zeros( (batch_size,self.pooled_chrom_window) )
        #chip_labels_batch = np.zeros( (batch_size,self.num_classes) )
        

        
        batch_start = epoch_tracker.cur_index
        batch_end = batch_start+batch_size
        batch_indices = self.perm_indices[batch_start:batch_end]

        #Defining decode_func here may help improve speed

        for bi, pindex in enumerate(batch_indices):
            
            (dna_seq_data,
             dna_shape_data,
             dnase_data,
             chip_labels_data ) = DreamDataSaverHdf5.decode_1d_to_nd(self.data[pindex],self.shape_list)
            
            
            self.dna_seq_batch[bi,0,:,:]=dna_seq_data
            self.dna_shape_batch[bi,0,:,:]=dna_shape_data
            self.dnase_seq_batch[bi,:] =dnase_data
            self.chip_labels_batch[bi,:] = chip_labels_data
            self.batch_indices [bi]= bi
        print "Pull/decode batch time is",time.clock()-t0
        #return (dna_seq_batch,dna_shape_batch,dnase_seq_batch,chip_labels_batch,batch_indices)


    def pull_batch_train(self,batch_size):
        #Pull a batch for training (calls training epoch tracker)
        #gc.disable()
        t0 = time.time()
        cProfile.runctx('self._pull_batch(batch_size,self.epoch_tracker)',globals(),locals())
        t1 = time.time()-t0
        print "_pull_batch time is",t1
        #Increment epoch tracker
        
        epoch_completed = self.epoch_tracker.increment(batch_size)
        
        if epoch_completed:
            
            #Re shuffle training indices if end of file reached.
            print "Starting on epoch",self.epoch_tracker.num_epochs
            self.perm_indices = np.random.permutation(self.perm_indices)

        #return self._pull_batch(batch_size,self.epoch_tracker)
        #return (self.dna_seq_batch[:],self.dna_shape_batch[:],self.dnase_seq_batch[:],self.chip_labels_batch[:],self.batch_indices[:])

    
    def pull_batch_eval(self,batch_size):
        # Pull a batch for evaluation. This calls the eval epoch tracker
        # And does not shuffle perm indices if end of file reached.
        # This is separate from pull_batch_train to avoid messing with
        # self.epoch_tracker.cur_index every time we need to pull items
        # for evaluation.
        dna_seq_batch,dna_shape_batch,dnase_seq_batch,chip_labels_batch,batch_indices =self._pull_batch(batch_size,self.eval_epoch_tracker)
        #Pull a batch for evaluation (calls evaluation epoch tracker)
        epoch_completed = self.eval_epoch_tracker.increment(batch_size)
        if epoch_completed:
            print "Finished evaluating an epoch"
        return dna_seq_batch,dna_shape_batch,dnase_seq_batch,chip_labels_batch,batch_indices
    




class DreamDataLoader:
    def __init__(self,
                 peaks_index_file,
                 coords_file,
                 cell_types,
                 genome_file,
                 genome_chromosome_sizes,
                 dna_shape_dir,
                 dnase_seq_dir,
                 rna_seq_dir,
                 genome_annotation_file,
                 chip_seq_labels_file,
                 seq_len = 600,
                 chrom_window = 3600,
                 chrom_pool_size=8):
                    

        self.peaks_index_file = peaks_index_file
        #Identify bound and unbound DNase peaks and use these for training

        self.indices = self.load_indices(self.peaks_index_file)
        self.num_examples =len(self.indices)
        print "Loading  selected indices from", self.indices
        print "Num examples: ", self.num_examples



        
        self.coords_file = coords_file.rstrip('/').rstrip('\\')
        print "Genome file is", genome_file.rstrip('/').rstrip('\\')
        self.dna_data = DnaData(genome_file.rstrip('/').rstrip('\\'))
        
        #Note:ChipSeqLabels will lookup celltypes from chip seq file if
        #  cell_types=None
        self.chip_seq_data = ChipSeqLabels(chip_seq_labels_file.rstrip('/').rstrip('\\'),
                                           cell_types) #Labels
        
        
        self.num_classes = 2 #This value must be two for binary classification

        #Cell types 
        self.cell_types = self.chip_seq_data.cell_types
        #if len(self.cell_types) ==1:
        #    self.cell_types = [self.cell_types]
        
        self.num_cell_types = len(self.cell_types)
        
        print "Cell types (input):", self.cell_types

        
        
        self.contig_size_dict = dbt.contig_sizes_to_dict(genome_chromosome_sizes)
        
        print "Init DNA shape data"
        self.dna_shape_data = DnaShapeData(dna_shape_dir)

        print "Init DNAse I chromatin accessibility data"
        self.chrom_window = chrom_window
        self.chrom_pool_size = chrom_pool_size
        self.dnase_data = DnaseData(dnase_seq_dir.rstrip('/').rstrip('\\'),
                                    self.cell_types,
                                    self.chrom_pool_size)
        print "Init RNA-seq data"
        self.rna_seq_data = RnaSeqData(rna_seq_dir,self.cell_types,
                                       genome_annotation_file)
        print "Init genome annotation"
        self.genome_annotation = GenomeAnnotationData(genome_annotation_file)             
        
        self.seq_len = seq_len
        
        self.indexer = DuGenCoords.IndexCoords(self.coords_file,
                                                     200,
                                                     50)
        
        self.collection = DataCollection(self,self.coords_file,
                                         self.indices,
                                         self.indexer,
                                         self.cell_types)
        
       
    
    def print_shuffled_indices(self,out_fname):
        with open(out_fname,'w') as out_file:
            for line in self.shuffled_indices:
                out_file.write(''.join(line))

    def load_indices(self,indices_fname):
        print "Loading indices file", indices_fname
        culled_indices = None
        with h5py.File(indices_fname) as hf:
                culled_indices = hf['selected_indices'][:]
        return culled_indices
    
class BaseData:
    def read_lines_to_list(self,n):
        with open(self.filename) as f:
            lines = [next(f) for x in xrange(n)]
        return lines
    
    @staticmethod
    def check_filename_ext(fname):
        return os.path.splitext(fname)[-1].lower()

    @staticmethod
    def is_gzip(fname):
        file_ext = BaseData.check_filename_ext(fname)
        if file_ext == '.gz' or file_ext == '.gzip':
            return True
        else:
            return False



class DreamDataLoaderJson(DreamDataLoader):
    #Initialize a DreamDataLoader object from a Json file
    def __init__(self, json_file, peaks_index_file,mode='train'):
        #Mode can be 'train','test','eval'
        
        params = DuJsonParse.JsonParams(json_file)
        if mode == 'train':
            DreamDataLoader.__init__(self,
                                peaks_index_file = peaks_index_file,
                                coords_file = params.train_coords_file,
                                cell_types = params.train_cell_types,
                                genome_file = params.genome_file,
                                genome_chromosome_sizes = params.genome_chromosome_sizes_file,
                                dna_shape_dir=params.dna_shape_dir,
                                dnase_seq_dir=params.dnase_seq_dir,
                                rna_seq_dir = params.rna_seq_dir,
                                genome_annotation_file= params.genome_annotation_file,
                                chip_seq_labels_file = params.chip_seq_labels_file,
                                seq_len = params.seq_len,
                                chrom_window = params.chrom_window,
                                chrom_pool_size = params.chrom_pool_size
                               )
        elif mode == 'test':
                DreamDataLoader.__init__(self,
                                peaks_index_file = peaks_index_file,
                                coords_file = params.test_coords_file,
                                cell_types = params.test_cell_types,
                                genome_file = params.genome_file,
                                genome_chromosome_sizes=params.genome_chromosome_sizes_file,
                                dna_shape_dir=params.dna_shape_dir,
                                dnase_seq_dir=params.dnase_seq_dir,
                                rna_seq_dir = params.rna_seq_dir,
                                genome_annotation_file= params.genome_annotation_file,
                                chip_seq_labels_file = params.chip_seq_labels_file,
                                seq_len = params.seq_len,
                                chrom_window = params.chrom_window,
                                chrom_pool_size = params.chrom_pool_size
                               )
     
        elif mode == 'eval':
                DreamDataLoader.__init__(self,
                                peaks_index_file = peaks_index_file,
                                coords_file = params.eval_coords_file,
                                cell_types = params.eval_cell_types,
                                genome_file = params.genome_file,
                                genome_chromosome_sizes = params.genome_chromosome_sizes_file,
                                dna_shape_dir = params.dna_shape_dir,
                                dnase_seq_dir = params.dnase_seq_dir,
                                rna_seq_dir = params.rna_seq_dir,
                                genome_annotation_file = params.genome_annotation_file,
                                chip_seq_labels_file = params.chip_seq_labels_file,
                                seq_len = params.seq_len,
                                chrom_window = params.chrom_window,
                                chrom_pool_size = params.chrom_pool_size
                               )
        else:
            print "Incorrect mode specified. 'train','test', and 'eval'."
        self.params = params
        


                       
class DnaData(BaseData):
    def __init__(self, filename = ''):
        self.filename = filename

        if BaseData.is_gzip(self.filename):
            print "Need to decompress genome fasta file", self.filename
        else:
            #Note: I ran tests to see how the speed of using pysam
            # compared to that of biopython. Using an indexed fasta file is
            # 3000X faster and vastly more memory efficient! 
            self.genome_idx = pysam.FastaFile(self.filename)
                

    def extract_coords_to_4d_onehot(self,contig_header,coord_lower,coord_upper):
        return dbt.seq_to_4d_onehot(self.extract_coords(contig_header,coord_lower,coord_upper))

    def extract_coords(self,contig_header,coord_lower,coord_upper):
        return self.genome_idx.fetch(contig_header,coord_lower,coord_upper)

    
class ChipSeqLabels(BaseData):
  
    def __init__(self, filename, cell_types, ambiguous_val =[0,0]):
        input_file_ext = os.path.splitext(filename)[1]
        if input_file_ext == '.gz' or input_file_ext == '.bgz':
            self.filename = filename
            print "ChIP-seq training file is", self.filename
        else:
            print "Input chip-seq training file must be sorted bgzip file"
            print "Use cmd: sort -k1,1 -k2,2n $filename|bgzip -c> $filename.bgz"
            print " followed by: tabix -p vcf $filename.bgz"
            self.filename = None
        self.pysam_tabix = pysam.TabixFile(self.filename)
        self.num_lines = dbt.count_lines_file(self.filename)-1 #Subtract 1 for header
        if cell_types != None:
            self.cell_types = cell_types.strip().split(',') #This should be a list
        else:
            print "Using cell types found in",self.filename
            self.cell_types = ChipSeqLabels.determine_cell_types(self.filename)

        #TODO: make it possible to only load specified cell types
        #Right now all cell types are loaded in the order they are listed
        

        print "Cell types are:",self.cell_types 
        self.num_cell_types = len(self.cell_types)
        self.label_dict = {'B':[0,1], 'A':ambiguous_val, 'U':[1,0]}
        #Note: May want to select different label for ambiguous label 'A'

        #self.train_labels = self.load_chip_labels()


    @staticmethod
    def determine_cell_types(chip_seq_labels_file):
        #Check first line of labels file to populate cell cell_types list
        #as well as determine number of classes for classification
        ext = os.path.splitext(chip_seq_labels_file)[1]
        if  ext == '.gz' or ext == '.bgz':
            cf = gzip.open(chip_seq_labels_file,'r')
        elif ext == '.tsv':
            cf = open(chip_seq_labels_file,'r')
        else:
            cf = None
            print "Error! Can\'t read label file:", chip_seq_labels_file
                
        header = cf.readline()
        #The header entries from header #4 onwards should be cell types
        cell_types = header.strip().split('\t')[3:]
        
        cf.close()
        return cell_types
   
    def load_coords_tabix(self,contig,start,end,cell_type_index):
        
        #Since tabix retreives only overlapping
        # queries, we only need to query an entry around the start coord
        query_window = 25
        query_start = start-query_window
        query_end = start+query_window

        try:
            records = self.pysam_tabix.fetch(contig,query_start,query_end)
            #records = dutabix.query(self.filename, contig,query_start,query_end)
        except:        
            print "Error in load_coords_tabix"
            print (contig,query_start,query_end,cell_type_index)
        #    #Making an arbitrary record so program won't crash
        #    records = self.pysam_tabix.fetch(contig,start-5,end+5)

        counter =0
        while True:
            
            try:
                rec = records.next()
                rec=rec.strip().split()
                counter += 1
                if int(rec[1]) == start and int(rec[2])==end:
                    query = rec
            except StopIteration:
                if counter==0:
                    print "No records detected in query region"
                    print "Setting query to a default value"
                    query = ['chr10', '575000', '575200', 'U','U','U','U','U','U']
                break
            
        #I had to put in a default value for query to stop the program from crashing
        #when records has nothing inside!!!
        #query = ['chr10', '575000', '575200', 'U']   
        #for rec in records:
        #    rec = rec.strip().split()
        #    if int(rec[1]) == start and int(rec[2])==end:
        #        query = rec
        #    else:
        #        print'\n'
        #        print 'Record is', rec
        #        print contig,':',start,'-',end,',',cell_type_index
        #        print "Exact coordinate match not found!"
        #        print "Tabix coordinate loading problem!"
        #        query = rec

        
        return self.label_dict[query[3+cell_type_index]]
        
        
class ChipSeqRelaxed(BaseData):
    name_struct = ['ChIPseq','cell_type','tf','relaxed','narrowPeak','sorted','bgz']
    def __init__(self,relaxed_peaks_dir,cell_types):
        #IDR score is column 6
        pass
        
    def load_coords_tabix(self,contig,start,end,cell_type_index):
        query_window = 25
        query_start = start-query_window
        query_end = start+query_window
        records = self.pysam_tabix.fetch(contig,query_start,query_end)
        for rec in records:
            pass
        
    

class DnaseData(BaseData):
    bw_file = ['DNASE','','fc','signal','bigwig']
    def __init__(self, directory,cell_types,chrom_pool_size):

        self.directory = directory
        self.chrom_pool_size = chrom_pool_size
        self.cell_types = cell_types
        self.filenames = self.get_filenames(cell_types)


        print "DNase I accessibility training files:", self.filenames
        


        
    def get_filenames(self,cell_types_list):
        filenames=[]
        for file in os.listdir(self.directory):
            cur_cell_file = os.path.basename(file).split('.')[1]
            if (cur_cell_file in cell_types_list):
                filenames.append(file)
        return filenames

                
    def load_single_coords(self,contig,start,end,cell_type,type='mean',feature_norm = 'standardize'):
        #print "In puller"
        #print (contig,start,end)
        #if (len(cell_type)==1):
        #    cell_type = [cell_type]
        fname_base = '.'.join([DnaseData.bw_file[0]]+[cell_type]+DnaseData.bw_file[2:])
        fname = self.directory+os.sep+fname_base
        bw = pyBigWig.open(fname)
        nBins = int((end-start)//self.chrom_pool_size)
        #Note: using bw.stats is faster than bw.values if you want to use bins
        # (I tested this myself!)
        vals= bw.stats(contig,start,end,type=type,nBins=nBins)
        if feature_norm=='minmax_rescale':
            #Perform min-max feature normalization
            vals = np.nan_to_num(vals)
            max_val = float(bw.header()['maxVal'])
            min_val= float(bw.header()['minVal'])
            try:
                return_vals = (vals-min_val)/(max_val-min_val)
            except TypeError:
                #print vals,max_val
                #print "vals has error in dnase pull"
                print "Attempting to convert None to 0 for dnase data"
                max_val = float(bw.header()['maxVal'])
                min_val= float(bw.header()['minVal'])
                med_val = 0.5
                vals = np.asarray([med_val if x==None else x for x in vals])
                
                return_vals = (vals-min_val)/(max_val-min_val)
                
        elif feature_norm == 'standardize':
            
            
            vals = np.nan_to_num(vals)
            #print "dnase sum data", bw.header()['sumData']
            #print "dnase n bases", bw.header()['nBasesCovered']
            mean_val = bw.header()['sumData']/float(bw.header()['nBasesCovered'])
            stddev  = np.sqrt(bw.header()['sumSquared']/float(bw.header()['nBasesCovered']))
            #print 'dnas stddev:',stddev
            #print 'dnase mean_val:',mean_val
                       
            try:
                return_vals = (vals-mean_val)/stddev
            except TypeError:
                #print vals, mean_val, stddev
                #print "vals has error in dnase pull"
                #print vals
                print "Dnase vals None detected. Attempting to convert None to mean"
                mean_val = bw.header()['sumData']/bw.header()['nBasesCovered']
                vals = np.asarray([mean_val if x==None else x for x in vals])
                stddev  = np.sqrt(bw.header()['sumSquared']/float(bw.header()['nBasesCovered']))
                return_vals = (vals-mean_val)/stddev
           
                
        else:
            return_vals= np.nan_to_num(
                          bw.stats(contig,start,end,type=type,nBins=nBins)
                                  )


        #return_val=bw.values(contig,start,end)
        
        bw.close()
        return return_vals

        
class DnaShapeData(BaseData):
    extensions=['.MGW','.Roll','ProT','HelT']
    extensions_bw = ['.MGW.wig.bw','.Roll.wig.bw','.ProT.wig.bw','.HelT.wig.bw']
   
    def __init__(self,dna_shape_dir):
        self.num_shapes = len(DnaShapeData.extensions)
        self.directory = dna_shape_dir
        
    @classmethod
    def init_single_file(cls,filename):
        self.num_shapes = len(DnaShapeData.extensions)
        self.filenames =[filename]
        
    
    def load_single_coords_bw(self, contig,start,end,feature_norm='standardize'):
        # Query dna dnase data for all filetypes listed in self.extensions
        # For four filetype (ie: MGW,Roll, ProT, and HelT)
        #  total return len will by a 4x*query_len numpy array
        # with rows place in order of listing (MGW,Roll,ProT, HelT)
        nuc_len =end-start
        #total_len =nuc_len*len(self.extensions) #Len of total return array
        #Preallocate return array
        all_stats = np.zeros( ( self.num_shapes,nuc_len) )
        bw_dir = os.listdir(self.directory)
        
        for i,ext in enumerate(DnaShapeData.extensions_bw):
            #print "Extension",ext
            for fname in bw_dir:
                #print self.directory+os.sep+fname
                fname = self.directory+os.sep+fname
                fname_ext = '.'+'.'.join(fname.split('.')[1:])
                if ext == fname_ext:
                    bw = pyBigWig.open(fname)
                    #stats = np.nan_to_num(bw.values(contig,start,end))
                    if feature_norm=='minmax_rescale':
                        max_val = float(bw.header()['maxVal'])
                        min_val = float(bw.header()['minVal'])
                        vals = np.nan_to_num(bw.values(contig,start,end))
                        all_stats[i,:] = (vals-min_val)/(max_val-min_val)
                    elif feature_norm == 'standardize':
                        mean_val = bw.header()['sumData']/float(bw.header()['nBasesCovered'])
                        
                        stddev  = np.sqrt(bw.header()['sumSquared']/
                                          float(bw.header()['nBasesCovered']))
                        vals = np.nan_to_num(bw.values(contig,start,end))
                        #print 'fname extension',fname_ext
                        #print 'dna shape stddev:',fname_ext,stddev
                        #print 'dna shape mean_val:',fname_ext,mean_val
                        all_stats[i,:] = (vals-mean_val)/stddev
                    else:
                        all_stats[i,:] = np.nan_to_num(bw.values(contig,start,end))
                    bw.close()
                #else:
                #print 'Error,', fname,'with extension', ext,'does not have extension', fname_ext
                ##print "Error, DnaShape couldn't be retrieved", (contig,start,end)
        #print all_stats
        return all_stats.T

    


class DnaShapeDataHdf5(BaseData):
    def __init__(self, data_dir=''):
        self.directory = data_dir
        self.types = ['MGW','Roll','ProT','HelT']
        self.extensions= [('.'+type) for type in self.types]
        self.num_shape_params = len(self.extensions)
        self.filenames = []
        for fname in os.listdir(self.directory):
            if fname.endswith('.h5'):
                self.filenames.append(file)

    def load_all_shape_coords(self, contig,start,end):
        #Load coordinates for all file extension.
        #Each vector will be appended in order to form a
        # long 1D vector containing all pertinent data
    
        nuc_len = len(start-end)
        #total_len =nuc_len*len(self.extensions) #Len of total return array
        #Preallocate return array
        all_stats = np.zeros((len(self.extensions),nuc_len))
    
        for ext in self.extensions:
            for fname in filenames:
                fext =os.path.splitext(os.path.basename(fname))[1]
                if fext == ext:
                    stat = self.load_single_coord(ext,contig,start,end)
                    all_stats[i,:] = stat
        return all_stats
                    
                
    def load_single_shape_coord(self,ext,contig,start,end):
        with h5py.File(fname,'r') as hf:
            return hf[contig][start:end]
                   
        
class RnaSeqData:
    def __init__(self, rna_seq_dir, cell_types,annotation_fname):
        self.directory = rna_seq_dir
        self.cell_types = cell_types

        self.annotation_filename = annotation_fname 

class GenomeAnnotationData:
    def __init__(self,genome_annotation_file):
        self.filename = genome_annotation_file
        
class DataCollection:
    def __init__(self, dream_loader, coords_file, indices, indexer, cell_types):

        self.dream_loader = dream_loader
        print "Cell types in collection", cell_types
        self.cell_types = cell_types
        self.num_cell_types = len(self.cell_types)
        #print "Num cell types", self.num_cell_types

        #Indices passed to constructor should already account for number
        # of cell types
        self.indices = np.array(indices)

        #The indexer generates coordinate data from a merged.bed coords_file
        # from, and can pull coordinate data using only an input index
        self.indexer = indexer
        

        self.num_examples= len(self.indices)
        #Examples is num_coords*num_cell_types, because data from
        #independent cell types is treated uniquely

        self.epoch_tracker = EpochTracker(self.num_examples)
        self.eval_epoch_tracker = EpochTracker(self.num_examples)
        self.coords_file = coords_file
        
        self.num_dna_shapes = 4
        
        #shuffle indices.
        self.perm_indices = np.random.permutation(self.indices) 
        '''
        Note: Always use perm_indices for training. perm_indices should
        be reshuffled with the completion of every epoch. 
        '''
        
        #Note: self.perm_indices will always reference the index of data within coord_data,
        # not the line index of data within the original input file

        #Note: For training instances indices should be shuffled prior
        # to training. This should be done before calling the
        # DataCollection constructor. Test/eval data does not necessarily 
        # need to be shuffled.        

        
        
    """
    ###DEPRECATED:     
    def load_coord_data(self,coords_file, num_cell_types,indices):

    # Retrieve coordinates with specified indices from coords file
    # if there are multiple cell types, repeat the loading for each cell type
    
        if (os.path.splitext(coords_file)[1] =='.gz'):
            f = gzip.open(coords_file,'r')
        elif (os.path.splitext(coords_file)[1] == '.bed'):
            f = open(coords_file,'r')
            
    
        #I nested this function for the sake of encapsulation
        num_lines = dbt.count_lines_file(coords_file)
        num_output_coords = len(indices)
        data_list = num_output_coords*num_cell_types*[None]#My attempt at preallocating
        for cell_ind in range(num_cell_types):
            for line
            _num,line in enumerate(f):
                if line_num in indices: #save only marked lines
                    contig,start,end = line.strip().split('\t')
                    start = int(start)
                    end = int(end)
                    data_list[cell_ind*num_output_coords+line_num] = (contig,
                                                                      start,
                                                                      end,
                                                                      cell_ind)
        f.close()
        return data_list
    """
    
    def _pull_batch(self,batch_size,epoch_tracker,pull_labels=True):
        dna_window = self.dream_loader.seq_len #600
        chrom_window = self.dream_loader.chrom_window #4800
        
        base_len = 200
        window_delta = int((dna_window-base_len)//2)
        chrom_window_delta = int((chrom_window-base_len)//(2))
        #window_delta describes how much to subtract or add to base coords
        #to get nucleotide window of nuc_window size
             
        #Preallocate return arrays
        #Treat different nucleotide letters as different channels
        dna_batch = np.zeros((batch_size,1,dna_window,4))

        
        dna_shape_batch = np.zeros((batch_size,
                                    1,
                                    dna_window,                                    
                                    self.dream_loader.dna_shape_data.num_shapes))


        ##Includes both dna sequence and dna shape 
        #dna_batch = np.zeros(batch_size,
        #                     1,
        #                     dna_window,
        #                     4+self.dream_loader.dna_shape_data.num_shapes)

    #Note: In dnase_seq_batch dim[1] is window size after average pooling
        dnase_seq_batch = np.zeros((batch_size,
                        int(chrom_window//self.dream_loader.dnase_data.chrom_pool_size) ))
        #rna_seq_batch = np.zeros(())
        #annotation_batch = np.zeros(())
        chip_labels_batch = np.zeros((batch_size,self.dream_loader.num_classes))
        
        batch_start = epoch_tracker.cur_index
        batch_end = batch_start+batch_size
        batch_indices = self.perm_indices[batch_start:batch_end]
        
        
        for bi,pindex in enumerate(batch_indices):
            #pindex is the shuffled index of the entry
            #bi is the index of the current element in the current batch

            contig,start,end,cell_index =  self.indexer.retrieve_by_index(pindex)   
           
            
            window_start = start-window_delta  #User set window. May be wider than 200 bp
            window_end   = end+window_delta      #User set window. May be wider than 200 bp
            #dna_batch shape should be [batch_ind,height=1,seq_len=dna_window,channels=4]
            #should concatenate on axis 0

            chrom_window_start = start-chrom_window_delta
            chrom_window_end   =  end+chrom_window_delta
            

            #If outside of chromsome bounds, move query window back within chromosome bounds
            contig_end = int(self.dream_loader.contig_size_dict[contig])

            if (chrom_window_start < 0):
                chrom_window_start = 0
                chrom_window_end = chrom_window_start+chrom_window
                window_start= int((chrom_window_start+chrom_window//2) - dna_window//2)
                window_end= window_start+dna_window
                print "chrom window queried start of chromosome"
            elif (chrom_window_end > contig_end):
                chrom_window_start = contig_end-chrom_window
                chrom_window_end = contig_end
                window_start= chrom_window_end - chrom_window//2 + dna_window//2
                window_end= window_start+dna_window
                
                print "chrom window queried end of chromosome"
            #elif (window_start<0) :
            #    window_start = 0
            #    window_end = window_start+dna_window
            #    print "dna window queried start of chromosome"
            #elif (window_end>contig_end):
            #    window_start = contig_end-dna_window
            #    window_end = contig_end
            #    print "dna window queried end of chromosome"
                
 #           t0 = time.clock()
            #DNA batch
            dna_batch[bi,0,:,:] = self.dream_loader.dna_data.extract_coords_to_4d_onehot(contig,
                                                              window_start,
                                                              window_end)
            #print "DNA batch retrieval time", time.clock()-t0
            #DNA shape batch
            dna_shape_batch[bi,0,:,:] = self.dream_loader.dna_shape_data.load_single_coords_bw(
                                                                      contig,
                                                                      window_start,
                                                                      window_end)
            #print "DNA shape batch retrieval time", time.clock()-t0
            #DNAse-seq batch
            #print ("contig",contig, "start",
            #       chrom_window_start, "end", chrom_window_end,"cell_ind",cell_index)
            dnase_seq_batch[bi,:] = self.dream_loader.dnase_data.load_single_coords(
                                       contig,
                                       chrom_window_start,
                                       chrom_window_end,
                                       self.cell_types[cell_index])
   #         print "DNAse seq batch retrieval time", time.clock()-t0
            #RNAse-seq batch
            #TODO

            #Genome annotation batch

            #ChIP labels with shape (batch_size,2)
            if pull_labels==True:
                chip_labels_batch[bi,:] = self.dream_loader.chip_seq_data.load_coords_tabix(contig,start,end,cell_index)
        #Apply relu to each element of label_set
            else:
                chip_labels_batch[bi,:]=[0.0,0.0]


           # print "ChIP batch retrieval time", time.clock()-t0
        return [dna_batch,dna_shape_batch, dnase_seq_batch,chip_labels_batch,batch_indices]



    def _pull_batch_b(self,batch_size):
        #Concatenate dna sequence and shape data into one matrix
        dna_batch,dna_shape_batch,dnase_seq_batch,chip_labels_batch,_ = self._pull_batch(batch_size)

        seq_shape_batch = np.concatenate((dna_batch,dna_shape_batch),axis=3)

        #print seq_shape_batch.shape
        return seq_shape_batch,dnase_seq_batch,chip_labels_batch
    
    
    def pull_batch_train(self,batch_size):
        #Pull a batch for training (calls training epoch tracker)        
        batch =self._pull_batch(batch_size,self.epoch_tracker)
            
        #Increment epoch tracker
        is_end_of_file = self.epoch_tracker.increment(batch_size)
        if is_end_of_file:
            #Re shuffle training indices if end of file reached. 
            self.perm_indices = np.random.permutation(self.perm_indices)
        return batch


    
    def pull_batch_eval(self,batch_size):
        # Pull a batch for evaluation. This calls the eval epoch tracker
        # And does not shuffle perm indices if end of file reached.
        # This is separate from pull_batch_train to avoid messing with
        # self.epoch_tracker.cur_index every time we need to pull items
        # for evaluation.
        batch =self._pull_batch(batch_size,self.eval_epoch_tracker)
        
        #Pull a batch for evaluation (calls evaluation epoch tracker)
        self.eval_epoch_tracker.increment(batch_size)
        return batch

    def pull_batch_prediction(self,batch_size):
        # Pull a batch for evaluation. This calls the eval epoch tracker
        # And does not shuffle perm indices if end of file reached.
        # This is separate from pull_batch_train to avoid messing with
        # self.epoch_tracker.cur_index every time we need to pull items
        # for evaluation.
        batch =self._pull_batch(batch_size,self.eval_epoch_tracker,pull_labels=False)
        
        #Pull a batch for evaluation (calls evaluation epoch tracker)
        self.eval_epoch_tracker.increment(batch_size)
        return batch


    
    def open(self):
        pass
    def close(self):
        pass
    



    
class EpochTracker:
    def __init__(self,num_examples):
        #Reminder: this exists as a seperate class to DataCollection
        #because the epoch tracking index need to be tracked separately during training
        # and evaluation
        self.num_examples = num_examples
        self.num_epochs = 0 #The number of epochs that have been passed
        self.cur_index = 0 #The index position on current epoch

    def increment(self,increment_size):
        #Returns true if end of current epoch
        new_index = self.cur_index + increment_size
        #Reset epoch counter if end of current epoch has been reached.
        if ( new_index >= self.num_examples):
            self.num_epochs += 1
            self.cur_index = 0
            #Reshuffle indices
            return True
        else:
            self.cur_index = new_index
            return False

        

          
      
    
