import numpy as np
import DuBioTools as dbt
import matplotlib.pyplot as plt
import os
from duseqlogo.LogoTools import LogoSheet

from duseqlogo.LogoTools import PwmTools

class TrainingReviewer():
    def __init__(self,
                sess,
                data_placeholder,
                labels_placeholder,
                keep_prob_placeholder,
                nuc_conv_model,
                nuc_data,
                save_dir='.'):
        self.session = sess
        self.data_placeholder = data_placeholder
        self.labels_placeholder = labels_placeholder
        self.keep_prob_placeholder = keep_prob_placeholder
        self.nuc_conv_model = nuc_conv_model
        self.nuc_data = nuc_data
        self.save_dir = save_dir

    def extract_simple_pwms(self,save_file = 'output.png',edge_clip=120):
        data = self.nuc_data.all_data
                                     
        #print 'Total data shape: ',data.shape
        h_conv_list = self._access_data(self.nuc_conv_model.h_conv1)

        #Retrieve max value index in each column (each filter)
        combined_pfm_list = self.nuc_conv_model.num_filters*[np.zeros((4,
                               self.nuc_conv_model.conv_filter_width),dtype=np.int)]
        #Visit every input data 
        for i,resp in enumerate(h_conv_list):

            #Examine each filter's response 
            for filt in range(resp.shape[3]): #resp.shape[3]=num_filters
                #Extract filter response, and squeeze dimensions that =1        
                filt_resp = np.squeeze(resp[:,:,:,filt])
                #filt_resp shape is [~seq_len]

                #plt.figure()
                #plt.plot(range(seq_len),cur_filter_resp)
                #plt.savefig(save_dir+os.sep+'test_'+str(j)+'_plot_test.png')

                #Extract subsequence with highest response to filter
                max_resp_ind = np.argmax(filt_resp)
                lower = int(1+max_resp_ind-self.nuc_conv_model.conv_filter_width//2) 
                upper = int(1+max_resp_ind+self.nuc_conv_model.conv_filter_width//2) 
                if lower>edge_clip and upper < (self.nuc_data.seq_len-edge_clip):
                    max_motif_onehot= np.squeeze(data[i,:,lower:upper,:],axis=0)
                    #print 'maxmotifshape', max_motif_onehot.shape
                    if max_motif_onehot.T.shape == combined_pfm_list[filt].shape:
                        combined_pfm_list[filt] = combined_pfm_list[filt] + max_motif_onehot.T.astype(int)
                    #Note for someone reason using the += operator in the above operation
                    #will cause a failure. I do not understand why.  
                    
                    #if filt==5:
                    #    data_slice = data[i,:,lower:upper,:]
                    #    dbt.onehot_4d_to_nuc(np.expand_dims(data_slice,axis=0))
                    
           
        #Now draw PWMs for each filter
        #print 'combined pfm', len(combined_pfm_list)
        #print combined_pfm_list[1]
        #print combined_pfm_list[2]
           
        logo_sheet = LogoSheet(combined_pfm_list,option='pwm',is_onehot_4d='True')
        logo_sheet.draw_pwm()
        logo_sheet.write_to_png(save_file)
        #print PwmTools.pfm_to_ppm(combined_pfm_list[5])
                
        #testlogos


    def _access_data(self,graph_node):
        data = self.nuc_data.all_data
        labels = self.nuc_data.all_labels
        num_examples = data.shape[0]
        seq_len = data.shape[2]
        num_filters = self.nuc_conv_model.num_filters
        conv_filter_width = self.nuc_conv_model.conv_filter_width
        #
        output_list = []

        for i in range(num_examples):
            node_output = self.session.run(graph_node,
                                             feed_dict=
                                              {self.data_placeholder:[data[i]],
                                               self.labels_placeholder:[labels[i]],
                                              self.keep_prob_placeholder:1.0})

            output_list.append(node_output)
        return output_list

    
        
