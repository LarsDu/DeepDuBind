import numpy as np
import DuBioTools as dbt

class NucDataCollection:
    #This class is used for hold-out cross validation
    def __init__(self,training_files='',test_files='',validation_files=''):
        #Convert to list if not a list
        self.training_files = self._check_if_list(training_files)
        self.test_files = self._check_if_list(test_files)
        self.validation_files= self._check_if_list(validation_files)
        #Print
        print 'Training files:\t',self.training_files
        print 'Test files:\t',self.test_files
        print 'Validation files:\t',self.validation_files
        #Extract data
        self.training_set = NucDataSet(self.training_files)
        self.test_set = NucDataSet(self.test_files)
        self.validation_set = NucDataSet(self.validation_files)

    def _check_if_list(self,filenames):
        if type(filenames) is not list:
            return [filenames]    
        else:
            return filenames


class NucKfoldCollection:
    """
    For k-fold cross validation, the starting NucDataSet needs to be
    is broken up into non-overlapping fractions. The hold_out fraction
    is used to determine test error. The average of all test errors is 
    will be used to calculate actual cross-validation error.
    The relevant NucData is stored under
    self.training_set[i] and self.validation_set[i]
    with different 'i' indices representing different ways of slicing the data
     
    """
    def __init__(self,file_list, k=3, validation_fraction = .2, shuffle=True):

        self.training_set = []
        self.validation_set = []

        self.validation_fraction = validation_fraction
        self.file_list = file_list
        self.all_data,self.all_labels,self.num_classes = dbt.extract_n_classes_fasta(file_list)
        print "Input files: ",file_list
        self.seq_len = self.all_data.shape[2]
        self.num_examples = self.all_data.shape[0]
        
        if (shuffle == True):
            #Shuffles self.data and self.labels
            self.perm = self.shuffle_all_data_labels() 
        else:
            self.perm = range(self.num_examples)
        self.validation_size = int(np.floor(self.validation_fraction*self.num_examples))
        self.training_size = int(self.num_examples-self.validation_size)
        print 'Training fraction size is: ', self.training_size,' examples'
        print 'Validation fraction size is: ', self.validation_size,' examples'
        
        self.validation_indices = []
        self.training_indices =[]

        #Populate training_set and validation_set lists
        #Each k index is a different way of slicing the data.
        for k_ind in range(k):
            self.validation_indices.append(
                self.perm[(self.validation_size*k_ind):(self.validation_size*(k_ind+1))]
                )
            #To get the remaining indices from self.perm, use np.setdiff1d
            self.training_indices.append(  np.setdiff1d(self.perm,
                                              self.validation_indices[k_ind]))
            
            
            cur_training_data = np.take(self.all_data,
                                        self.training_indices[k_ind],
                                        axis =0)
            #print "SHAPE!",cur_training_data.shape
            
            cur_training_labels = np.take(self.all_labels,
                                          self.training_indices[k_ind],
                                          axis =0)
            
            cur_validation_data = np.take(self.all_data,
                                          self.validation_indices[k_ind],
                                          axis=0)
            cur_validation_labels = np.take(self.all_labels,
                                            self.validation_indices[k_ind],
                                            axis=0)
            self.training_set.append(
                NucDataSet(cur_training_data,       
                        cur_training_labels,
                        self.num_classes,
                            shuffle=False)
                         )

            self.validation_set.append(
                NucDataSet(cur_validation_data,
                            cur_validation_labels,
                            self.num_classes,
                            shuffle=False)
                         )
            
                        
            
    def shuffle_all_data_labels(self):
        #Check that the number of examples match

        if self.all_data.shape[0] != self.all_labels.shape[0]:
            print 'Error, data and labels contain different number of examples'
            return None
        perm = np.arange(self.num_examples)
        np.random.shuffle(perm) #Shuffle perm in place
        self.all_data = np.take(self.all_data,perm,axis=0)
        self.all_labels = np.take(self.all_labels,perm,axis=0)
        return perm

    def write_summary_file():
        '''
        Write a file that summarizes
        '''
            


    

        
class NucDataSet:
    '''This class is for encapsulation of nucleotide data and labels
    in the form of one-hot vectors for input into tensorflow models.
    It also contains methods for shuffling the data and tracking
    epochs (an epoch is one training run involving an entire dataset).

    For long training runs, this object will reshuffle the data
    at the beginning of new epochs, while also keeping track of the
    number of epochs.

    Args:
    file_list: List of fasta files. Each file will get a one-hot label
    corresponding to order of input.
  
    '''
    # The data should be shuffled once. Once the entire data set has
    # been pulled (one epoch), it should be reshuffled for subsequent
    # batch_pulls
    def __init__(self,data,labels,num_classes, shuffle=True):
        
        self.file_list = ['']
        self.data = data
        self.labels = labels
        self.num_classes = num_classes
        self.num_examples = self.data.shape[0]
        
        if (shuffle==True):
            #Shuffle self.data and self.labels 
            self.perm = self.shuffle_data_labels() 
            #Counts the number of training examples pulled on the current epoch.
            #This counter should reset on each new epoch
        else:
            self.perm = range(self.num_examples)
        self.epoch_training_index = 0
        #Epoch counter
        self.num_epochs = 0

    @classmethod
    def init_from_files(cls,file_list):
        '''
        Alternate to default constructor
        Initialize self.data,self.labels,self.num_classes from a
        list of files
        '''
        self.file_list = file_list
        self.data,self.labels,self.num_classes = dbt.extract_n_classes_fasta(file_list)
        self.num_examples = self.data.shape[0]
        self.perm =[] #Records the indices used for initial round of shuffling.
        if (shuffle==True):
            #Shuffle self.data and self.labels 
            self.perm = self.shuffle_data_labels() 
            #Counts the number of training examples pulled on the current epoch.
            #This counter should reset on each new epoch
        self.epoch_training_index = 0
        #Epoch counter
        self.num_epochs = 0


    
    def pull_batch(self,batch_size):
        self.epoch_training_index += batch_size
        #Check that on the next iteration, the epoch index
        #does not exceed the number of examples.
        #If this is the case, reset all counters, and
        #increment the epoch counter             
         #Pull a batch of previously shuffled data and labels
        data_batch = self.data[self.epoch_training_index:
                (self.epoch_training_index+batch_size)]
        labels_batch = self.labels[self.epoch_training_index:
            (self.epoch_training_index+batch_size)]

        #Check if the next batch overflows, and reset counters if it does
        if (self.epoch_training_index+batch_size) >= self.num_examples:
            #If an epoch has been completed, reshuffle the data, reset
            #the epoch index counter, and increment the epoch counter
            self.epoch_training_index=0
            self.num_epochs += 1
            self.shuffle_data_labels() #Shuffle self.data and self.labels
       
        return data_batch,labels_batch


    def shuffle_data_labels(self):
        #Check that the number of examples match
        if self.data.shape[0] != self.labels.shape[0]:
            print 'Error, data and labels contain different number of examples'
            return None
        perm = np.arange(self.num_examples)
        np.random.shuffle(perm) #Shuffle perm in place
        self.data = np.take(self.data,perm,axis=0)
        self.labels = np.take(self.labels,perm,axis=0)
        return perm

    def reset(self):
        self.shuffle_data_labels()
        self.epoch_training_index =0
        self.num_epochs=0
