import tensorflow as tf
import numpy as np
import time

class NucConvModel:

    def __init__(self,seq_len,mean_pooled_dnase_window,**kwargs):

        #These class variables are
        self.seq_len = seq_len #Recommend 600 bp
        self.mean_pooled_dnase_window = mean_pooled_dnase_window
        self.num_classes = 2
        self.dna_conv_filter_width = 24
        self.num_dna_filters = 48

        #Note the chromatin_window may be pooled on batch_pull, so
        # it does not reflect actual # bp's.
        # I've set it to be a 4800 bp window divided by 8, since resolution
        # is not a big deal here.


        #Number of dna shapes (ie: Roll, MGW, ProT, HelT)
        self.num_dna_shapes =4
        self.dna_pool_size = 4 #Maxpooling for conv layer
        self.num_fc1_neurons = 512
        self.num_fc2_neurons = 512
        self.num_fc3_neurons = 256
        self.num_fc4_neurons = 128
        self.num_fc5_neurons = 64
        
        #Specific params for inferenceA and inferenceC
        self.num_dnase_weights = 256
        self.dna_shape_multiplier = 8 #the num filters per dna shape (like MGW)
        self.dna_shape_conv_filter_width = 32

        #Params for inferenceC
        
        self.dna_shape_pool_size = 4
        self.dnase_pool_size =2
        self.dnase_conv_filter_width = 64
        self.num_dnase_conv_filters = 54

        
        #Option to overwrite defaults with keyword args
        #from constructor call

        #kwargs = sorted(kwargs.keys())
        for kw in kwargs:
            if kw == 'seq_len':
                self.seq_len = kwargs[kw]
                print "seq_len is now", self.seq_len
            elif kw == 'mean_pooled_dnase_window':
                self.mean_pooled_dnase_window = kwargs[kw]
                print "Mean pooled chromatin_window is:", self.mean_pooled_dnase_window
            elif kw == 'num_classes':
                self.num_classes = kwargs[kw]
            elif kw == 'conv_filter_width' :
                self.conv_filter_width = kwargs[kw]
            elif kw == 'num_filters':
                self.num_dna_filters = kwargs[kw]
            elif kw == 'pool_size':
                self.dna_pool_size = kwargs[kw]
            elif kw == 'num_fc1_neurons':
                self.num_fc1_neurons = kwargs[kw]
                print "num_fc1_neurons:", self.num_fc1_neurons
            elif kw == 'num_fc2_neurons':
                self.num_fc2_neurons = kwargs[kw]
                print "num_fc2_neurons:", self.num_fc2_neurons
            elif kw == 'num_fc3_neurons':
                self.num_fc3_neurons = kwargs[kw]
                print "num_fc3_neurons:", self.num_fc3_neurons
            elif kw == 'num_fc4_neurons':
                self.num_fc4_neurons = kwargs[kw]
                print "num_fc4_neurons:", self.num_fc3_neurons
            elif kw == 'num_fc5_neurons':
                self.num_fc5_neurons = kwargs[kw]
                print "num_fc5_neurons:", self.num_fc3_neurons
        
                
            else:
                print 'No custom model parameters specified'

        #These class variables are for ease of data extraction,
        # and get computed in inference. They are listed
        # here for ease of reference
        self.dna_conv_filters = None #Calculated by inference
        self.fully_connected_weights1 = None #Set in inferenceB
        self.fully_connected_weights2 = None #Set in inferenceB
        self.fully_connected_weights3 = None #Set in inferenceB
        self.h_conv1 = None #Set in inference



    def inferenceA(self,dna_seq_data,dna_shape_data,dnase_data, keep_prob):
        ##Conv1
        with tf.variable_scope('dna_conv1') as scope:
            #Initialize convolution filters weights and output bias variables
            # Dims are [height,width,num_channels,num_filters]
            #Input is shape [batch_size,1,nuc_len,4] where dim[3]
            # represents different channels/nucleotide letters
            W_conv1 = init_weights('filter_weights',
                                   shape=[1,self.dna_conv_filter_width,4,self.num_dna_filters])

            self.dna_conv_filters = W_conv1
            b_conv1 = init_bias('biases',[self.num_dna_filters])

            # Apply convolution and ReLu to raw DNA data
            h_conv1 = tf.nn.relu(tf.nn.conv2d(input=dna_seq_data,filter=W_conv1,
                                            strides=[1,1,1,1],padding='SAME',
                                                name=scope.name)+b_conv1)
            self.h_conv1 = h_conv1
            activation_summary(h_conv1)

        ##MaxPool1    
        with tf.variable_scope('dna_max_pool1') as scope:
            h_pool1_max = tf.nn.max_pool(h_conv1, ksize = [1,1,self.dna_pool_size,1],
                                                strides = [1,1,self.dna_pool_size,1],
                                                padding = 'SAME',
                                                name = scope.name)

            h_pool1_max = tf.squeeze(h_pool1_max,[1])
            #h_poo1_max now has shape
            #[batch_size,self.seq_len/self.dna_pool_size, self.num_dna_filters]     
            pooled_elem_len = self.seq_len/self.dna_pool_size
            #Flatten the data into [batch_size, pooled_elem_len*self.num_dna_filters]
            ###FIX SQUEEZE
            h_pool1_max = tf.reshape(h_pool1_max,
                                            [-1,pooled_elem_len*self.num_dna_filters])


        #DNA shape convolution (no pooling!)
        with tf.variable_scope('dna_shape_conv1') as scope:
            #Convolve DNA shapes once

            #Weights example:[height=1,width=4,channels=4,multiplier=8]
            W_shape_conv1 = init_weights('dna_shape_filter_weights',
                                         shape=[1,
                                                self.dna_shape_conv_filter_width,
                                                self.num_dna_shapes,
                                                self.dna_shape_multiplier])
        
            b_shape_conv1 = init_bias('biases',[self.dna_shape_multiplier*self.num_dna_shapes])

            # Note: input should be shape [batch_size,1,nuc_len,4].
            # dim[3] holds MGW,Roll,ProT,HelT indices

            #Note: depthwise_conv2d applies a different filter for each input channel
            # In this case each input channel is a DNA shape (MGW,Roll,etc)
            # The results are concatenated together
            
            h_shape_conv1 = tf.nn.relu(
                              tf.nn.depthwise_conv2d(input=dna_shape_data,filter=W_shape_conv1, strides=[1,1,1,1], padding='SAME',name=scope.name)+b_shape_conv1)
            """h_shape_conv1 should end up with shape [batch_size, height,width,
               num_filters_per_channel*num_channels] which in this context should
               be [batch_size,1,nuc_len,dna_shape_multiplier*4]                   """
            h_shape_conv1 = tf.squeeze(h_shape_conv1,[1])
            #Flatten data to shape [batch_size, filtered_data]
            h_shape_conv1 = tf.reshape(h_shape_conv1,
                    [-1,self.seq_len*self.dna_shape_multiplier*self.num_dna_shapes])


            activation_summary(h_shape_conv1)

            
        
        with tf.variable_scope('dnase') as scope:
            W_dnase1 = init_weights('dnase_weights',
                                    shape=[self.mean_pooled_dnase_window,
                                            self.num_dnase_weights],
                                            mean=0)
            b_dnase1= init_bias('biases',[self.num_dnase_weights])
            h_dnase1 = tf.nn.relu(tf.matmul(dnase_data,W_dnase1)+b_dnase1)
            
        ##Fully-connected layer 1
        with tf.variable_scope('fc1') as scope:
            #Determine total length of 1d vector being connected to fully
            #connected layers 
            dna_elem_len = pooled_elem_len*self.num_dna_filters
            dna_shape_elem_len = self.seq_len*self.dna_shape_multiplier*self.num_dna_shapes
            dnase_elem_len = self.num_dnase_weights
            total_elem_len = dna_elem_len+dna_shape_elem_len+dnase_elem_len


            
            W_fc1= init_weights('weights',[total_elem_len,
                                    self.num_fc1_neurons])
            h_fc1 = tf.concat(1,[h_pool1_max,h_shape_conv1,h_dnase1])
            self.fully_connected_weights1  = W_fc1
            b_fc1 = init_bias('biases',[self.num_fc1_neurons])
            h_fc1 = tf.nn.relu(tf.matmul(h_fc1,W_fc1)+b_fc1, name=scope.name)
            activation_summary(h_fc1)

        #Fully connected layer 2
        with tf.variable_scope('fc2') as scope:    
            W_fc2 = init_weights('weights',[self.num_fc1_neurons,
                                            self.num_fc2_neurons])
            self.fully_connected_weights2 = W_fc2
            b_fc2 = init_bias('biases',[self.num_fc2_neurons])
            h_fc2 = tf.nn.relu(tf.matmul(h_fc1,W_fc2)+b_fc2,name=scope.name)
            activation_summary(h_fc2)

        ## Dropout - Apply dropout before readout layer to reduce overfitting
        with tf.variable_scope('dropout') as scope:
            #keep_prob = tf.placeholder(tf.float32,name="keep_prob")
            h_fc2_drop = tf.nn.dropout(h_fc2, keep_prob,name=scope.name)

        # Readout
        with tf.variable_scope('readout_logits') as scope:
            W_fcr = init_weights('weights',[self.num_fc2_neurons,self.num_classes])
            b_fcr = init_bias('biases',[self.num_classes])
            logits=tf.matmul(h_fc2_drop, W_fcr) + b_fcr
            activation_summary(logits)
        return logits



    
    def inferenceB(self,dna_seq_data,dna_shape,dnase_data,keep_prob):
        ##Conv1

        seq_shape_data = tf.concat(3,[dna_seq_data,dna_shape])
        with tf.variable_scope('seq_shape_conv1') as scope:
            #Initialize convolution filters weights and output bias variables
            # Dims are [height,width,num_channels,num_filters]
            #Input is shape [batch_size,1,nuc_len,4+num_shapes]
            # dim[3] is the four nucleotide letter followed by the four
            # dna shape parameters (MGW,Roll, HelT, ProT)
            W_conv1 = init_weights('filter_weights',
                                   shape=[1,
                                          self.dna_conv_filter_width,
                                          4+self.num_dna_shapes,
                                          self.num_dna_filters])

            
            b_conv1 = init_bias('biases',[self.num_dna_filters])

            # Apply convolution and ReLu to raw DNA data
            
            h_conv1 = tf.nn.relu(tf.nn.conv2d(input=seq_shape_data,filter=W_conv1,
                                            strides=[1,1,1,1],padding='SAME',
                                                name=scope.name)+b_conv1)
            self.dna_conv_filters = W_conv1
            self.h_conv1 = h_conv1
            activation_summary(h_conv1)

        ##MaxPool1    
        with tf.variable_scope('seq_shape_max_pool1') as scope:
            h_pool1_max = tf.nn.max_pool(h_conv1, ksize = [1,1,self.dna_pool_size,1],
                                                strides = [1,1,self.dna_pool_size,1],
                                                padding = 'SAME',
                                                name = scope.name)

            h_pool1_max = tf.squeeze(h_pool1_max,[1])
            #h_poo1_max now has shape
            #[batch_size,self.seq_len/self.dna_pool_size, self.num_dna_filters]     
            pooled_elem_len = self.seq_len/self.dna_pool_size
            #Flatten the data into [batch_size, pooled_elem_len*self.num_dna_filters]

            ###FIXSQUEEZE
            h_pool1_max = tf.reshape(h_pool1_max,
                                            [-1,pooled_elem_len*self.num_dna_filters])

        with tf.variable_scope('fc1') as scope:
            W_dna = init_weights('weights',
                        [pooled_elem_len*self.num_dna_filters,
                         self.num_fc1_neurons])

            #The mean for the dnase weights may be 0.5 instead of 0 so I've separated weight initialization here
            W_dnase = W_fc1 = init_weights('weights',
                        [self.mean_pooled_dnase_window,
                         self.num_fc1_neurons],mean=0)
            W_fc1 = tf.concat(0,[W_dna,W_dnase])
            #DNAse data injected here
            x_fc1 = tf.concat(1,[h_pool1_max,dnase_data])
            b_fc1 = init_bias('biases',[self.num_fc1_neurons])
            h_fc1 = tf.nn.relu(tf.matmul(x_fc1,W_fc1)+b_fc1,name=scope.name)
            
            self.fully_connected_weights1 = W_fc1
            activation_summary(h_fc1)

        with tf.variable_scope('fc2') as scope:
            W_fc2 = init_weights('weights',
                                 [self.num_fc1_neurons,self.num_fc2_neurons])
            b_fc2 = init_bias('biases',[self.num_fc2_neurons])
            h_fc2 = tf.nn.relu(tf.matmul(h_fc1,W_fc2)+b_fc2, name=scope.name)

            self.fully_connected_weights2 = W_fc2
            activation_summary(h_fc2)
            
        with tf.variable_scope('fc3') as scope:
            W_fc3 = init_weights('weights',
                                 [self.num_fc2_neurons,self.num_fc3_neurons])
            b_fc3 = init_bias('biases',[self.num_fc3_neurons])
            h_fc3 = tf.nn.relu(tf.matmul(h_fc2,W_fc3)+b_fc3,name=scope.name)

            self.fully_connected_weights3 = W_fc3
            activation_summary(h_fc3)

        with tf.variable_scope('fc4') as scope:
            W_fc4 = init_weights('weights',
                                 [self.num_fc3_neurons,self.num_fc4_neurons])
            b_fc4 = init_bias('biases',[self.num_fc4_neurons])
            h_fc4 = tf.nn.relu(tf.matmul(h_fc3,W_fc4)+b_fc4,name=scope.name)

            self.fully_connected_weights4 = W_fc4
            activation_summary(h_fc4)


        ## Dropout - Apply dropout before readout layer to reduce overfitting
        with tf.variable_scope('dropout') as scope:
            #keep_prob = tf.placeholder(tf.float32,name="keep_prob")
            h_fc4_drop = tf.nn.dropout(h_fc4, keep_prob,name=scope.name)

            
        #Readout
        with tf.variable_scope('readout_logits') as scope:
            W_fcr = init_weights('weights',[self.num_fc4_neurons,self.num_classes])
            b_fcr = init_bias('biases',[self.num_classes])
            logits=tf.matmul(h_fc4_drop, W_fcr) + b_fcr
            activation_summary(logits)
        return logits



    def inferenceC(self,dna_seq_data,dna_shape_data,dnase_data, keep_prob):
        """Convolution is peformed separately for dna_seq,dna_shape,
        and dnase data:

            - For dna sequence, each nucleotide occupies a channel, and
              the results of convolution passes are summed up across channels.
              A single maxpooling operation of 4 bp is then performed
              
            - For dna shape, separate convolutions are performed for each shape,
              and the convolved outputs are not added together.
              A single maxpooling operation of 4 bp is then performed on each output 
              
            - For dnase_data, a large linear convolution encompasing ~510 bp
              (85 post pooled ) is used, followed by a max pooling for every 4 units 
              
              The output of the three aforementioned operations is then flattened into
              a single vector and passed to a 4 layer fully connected network consisting
              of 512,512,256,64 neurons

        """

        
        ##Conv1
        with tf.variable_scope('dna_conv1') as scope:
            #Initialize convolution filters weights and output bias variables
            # Dims are [height,width,num_channels,num_filters]
            #Input is shape [batch_size,1,nuc_len,4] where dim[3]
            # represents different channels/nucleotide letters
            W_conv1 = init_weights('filter_weights',
                                   shape=[1,self.dna_conv_filter_width,4,self.num_dna_filters])

            self.dna_conv_filters = W_conv1
            b_conv1 = init_bias('biases',[self.num_dna_filters])

            # Apply convolution and ReLu to raw DNA data
            h_conv1 = tf.nn.relu(tf.nn.conv2d(input=dna_seq_data,filter=W_conv1,
                                            strides=[1,1,1,1],padding='SAME',
                                                name=scope.name)+b_conv1)
            self.h_conv1 = h_conv1
            activation_summary(h_conv1)

        ##DNA MaxPool1 (inference C)   
        with tf.variable_scope('dna_max_pool1') as scope:
            h_pool1_max = tf.nn.max_pool(h_conv1, ksize = [1,1,self.dna_pool_size,1],
                                                strides = [1,1,self.dna_pool_size,1],
                                                padding = 'SAME',
                                                name = scope.name)

            h_pool1_max = tf.squeeze(h_pool1_max,[1])
            #h_poo1_max now has shape
            #[batch_size,self.seq_len/self.dna_pool_size, self.num_dna_filters]     
            pooled_elem_len = self.seq_len/self.dna_pool_size
            #Flatten the data into [batch_size, pooled_elem_len*self.num_dna_filters]
            h_pool1_max_len = pooled_elem_len*self.num_dna_filters
            #Flatten
            ###FIX SQUEEZE
            h_pool1_max = tf.reshape(h_pool1_max,
                                            [-1,h_pool1_max_len])

            


        #DNA shape convolution (no pooling!) (inference C)
        with tf.variable_scope('dna_shape_conv1') as scope:
            #Convolve DNA shapes once

            #Weights example:[height=1,width=4,channels=4,multiplier=8]
            W_shape_conv1 = init_weights('dna_shape_filter_weights',
                                         shape=[1,
                                                self.dna_shape_conv_filter_width,
                                                self.num_dna_shapes,
                                                self.dna_shape_multiplier])
        
            b_shape_conv1 = init_bias('biases',[self.dna_shape_multiplier*self.num_dna_shapes])


            
            # Note: input should be shape [batch_size,1,nuc_len,4].
            # dim[3] holds MGW,Roll,ProT,HelT indices

            #Note: depthwise_conv2d applies a different filter for each input channel
            # In this case each input channel is a DNA shape (MGW,Roll,etc)
            # The results are concatenated together
            
            h_shape_conv1 = tf.nn.relu(
                              tf.nn.depthwise_conv2d(input=dna_shape_data,filter=W_shape_conv1, strides=[1,1,1,1], padding='SAME',name=scope.name)+b_shape_conv1)
            """h_shape_conv1 should end up with shape [batch_size, height,width,
               num_filters_per_channel*num_channels] which in this context should
               be [batch_size,1,seq_len,dna_shape_multiplier*4]                   """


        #DNA shape maxpooling (inference C)
        with tf.variable_scope('dna_shape_max_pool1') as scope:
            h_shape_max1 = tf.nn.max_pool(h_shape_conv1,
                                          ksize=[1,1,self.dna_shape_pool_size,1],
                                          strides=[1,1,self.dna_shape_pool_size,1],
                                          padding='SAME',
                                          name = scope.name)

            #Remove dim[1] which is used for height in 2d convolutions
            h_shape_max1 = tf.squeeze(h_shape_max1,[1])
            """
            h_shape_max1 should go from shape (batch_size,1,seq_len/dna_shape_pool_size,dna_shape_multiplier*4)
            to (batch_size,1,seq_len/dna_shape_pool_size,dna_shape_multiplier*4)
            """
            
            h_shape_max1_len = int((self.seq_len/self.dna_shape_pool_size)
                                   *self.dna_shape_multiplier
                                   *self.num_dna_shapes)

            h_shape_max1 = tf.reshape(h_shape_max1,[-1,h_shape_max1_len])
            activation_summary(h_shape_conv1)

            
        # DNase-Seq convolution (inferenceC)
        with tf.variable_scope('dnase_conv1') as scope:
            #Reshape data such that shape is now (batch_size,1,mean_pooled_dnase_window,1)
            dnase_data = tf.expand_dims(dnase_data,2)
            dnase_data = tf.expand_dims(dnase_data,1)
            
            


            W_dnase_conv1 = init_weights('dnase_weights',
                                    shape=[1,self.dnase_conv_filter_width,
                                           1,self.num_dnase_conv_filters])
            b_dnase_conv1= init_bias('biases',[self.num_dnase_conv_filters])
            h_dnase_conv1 = tf.nn.relu(tf.nn.conv2d(input=dnase_data,
                                              filter=W_dnase_conv1,
                                              strides=[1,1,1,1],
                                              padding='SAME',
                                              name=scope.name)+b_dnase_conv1)

        with tf.variable_scope('dnase_max_pool1') as scope:
            """Maxpool h_dnase_conv1 from shape (batch_size,1,mean_pooled_dnase_window,self.num_dnase_filters)
               to shape ((batch_size,1,mean_pooled_dnase_window/dnase_pool_size,self.num_dnase_filters)
            """
            
            h_dnase_max1 = tf.nn.max_pool(h_dnase_conv1,
                                         ksize = [1,1,self.dnase_pool_size,1],
                                         strides=[1,1,self.dnase_pool_size,1],
                                         padding='SAME',
                                         name=scope.name)
            
            h_dnase_max1=tf.squeeze(h_dnase_max1,[1])
            h_dnase_max1_len = int((self.mean_pooled_dnase_window/self.dnase_pool_size)*
                                   self.num_dnase_conv_filters)
            h_dnase_max1 = tf.reshape(h_dnase_max1,[-1,h_dnase_max1_len])

            
        ##Fully-connected layer 1
        with tf.variable_scope('fc1') as scope:
            #Determine total length of 1d vector being connected to fully
            #connected layers 
            total_elem_len = h_pool1_max_len+h_shape_max1_len+h_dnase_max1_len
            print "Total fc1 input vector length:", total_elem_len

            print h_pool1_max
            print h_shape_max1
            print h_dnase_max1
            W_fc1= init_weights('weights',[total_elem_len,
                                    self.num_fc1_neurons])
            h_fc1 = tf.concat(1,[h_pool1_max,h_shape_max1,h_dnase_max1])
            self.fully_connected_weights1  = W_fc1
            b_fc1 = init_bias('biases',[self.num_fc1_neurons])
            h_fc1 = tf.nn.relu(tf.matmul(h_fc1,W_fc1)+b_fc1, name=scope.name)
            activation_summary(h_fc1)

        #Fully connected layer 2
        with tf.variable_scope('fc2') as scope:    
            W_fc2 = init_weights('weights',[self.num_fc1_neurons,
                                            self.num_fc2_neurons])
            self.fully_connected_weights2 = W_fc2
            b_fc2 = init_bias('biases',[self.num_fc2_neurons])
            h_fc2 = tf.nn.relu(tf.matmul(h_fc1,W_fc2)+b_fc2,name=scope.name)
            activation_summary(h_fc2)

        with tf.variable_scope('fc3') as scope:
            W_fc3 = init_weights('weights',
                                 [self.num_fc2_neurons,self.num_fc3_neurons])
            b_fc3 = init_bias('biases',[self.num_fc3_neurons])
            h_fc3 = tf.nn.relu(tf.matmul(h_fc2,W_fc3)+b_fc3,name=scope.name)

            self.fully_connected_weights3 = W_fc3
            activation_summary(h_fc3)

        with tf.variable_scope('fc4') as scope:
            W_fc4 = init_weights('weights',
                                 [self.num_fc3_neurons,self.num_fc4_neurons])
            b_fc4 = init_bias('biases',[self.num_fc4_neurons])
            h_fc4 = tf.nn.relu(tf.matmul(h_fc3,W_fc4)+b_fc4,name=scope.name)

            self.fully_connected_weights4 = W_fc4
            activation_summary(h_fc4)


        ## Dropout - Apply dropout before readout layer to reduce overfitting
        with tf.variable_scope('dropout') as scope:
            #keep_prob = tf.placeholder(tf.float32,name="keep_prob")
            h_fc4_drop = tf.nn.dropout(h_fc4, keep_prob,name=scope.name)

            
        #Readout
        with tf.variable_scope('readout_logits') as scope:
            W_fcr = init_weights('weights',[self.num_fc4_neurons,self.num_classes])
            b_fcr = init_bias('biases',[self.num_classes])
            logits=tf.matmul(h_fc4_drop, W_fcr) + b_fcr
            activation_summary(logits)
        return logits





        
            

    def inference_nuc_only(self,nuc_data,keep_prob):
        ##Conv1
        with tf.variable_scope('conv1') as scope:
            #Initialize convolution filters weights and output bias variables
            # Dims are [height,width,num_channels,num_filters]

            W_conv1 = init_weights('filter_weights',
                                   shape=[1,self.dna_conv_filter_width,4,self.num_dna_filters])

            self.conv_filters = W_conv1
            b_conv1 = init_bias('biases',[self.num_dna_filters])

            # Apply convolution and ReLu
            h_conv1 = tf.nn.relu(tf.nn.conv2d(input=nuc_data,filter=W_conv1,
                                            strides=[1,1,1,1],padding='SAME',
                                                name=scope.name)+b_conv1)
            self.h_conv1 = h_conv1
            activation_summary(h_conv1)

        ##MaxPool1    
        with tf.variable_scope('dna_max_pool1') as scope:
            h_pool1_max = tf.nn.max_pool(h_conv1, ksize = [1,1,self.dna_pool_size,1],
                                                strides = [1,1,self.dna_pool_size,1],
                                                padding = 'SAME',
                                                name = scope.name)

            h_pool1_max = tf.squeeze(h_pool1_max,[1])
            #h_poo1_max now has shape
            #[batch_size,self.seq_len/self.dna_pool_size, self.num_dna_filters]     
            pooled_elem_len = self.seq_len/self.dna_pool_size 
            #Flatten the data into [batch_size, pooled_elem_len*self.num_dna_filters]
            ###FIX SQUEEZE
            h_pool1_max = tf.reshape(h_pool1_max,
                                            [-1,pooled_elem_len*self.num_dna_filters])
            
            ##Fully-connected layer 1
        with tf.variable_scope('fc1') as scope:    
            W_fc1= init_weights('weights',[pooled_elem_len*self.num_dna_filters,
                                    self.num_fc1_neurons])
            self.fully_connected_weights1  = W_fc1
            b_fc1 = init_bias('biases',[self.num_fc1_neurons])
            h_fc1 = tf.nn.relu(tf.matmul(h_pool1_max,W_fc1)+b_fc1, name=scope.name)
            activation_summary(h_fc1)

        #Fully connected layer 2
        with tf.variable_scope('fc2') as scope:    
            W_fc2 = init_weights('weights',[self.num_fc1_neurons,
                                            self.num_fc2_neurons])
            self.fully_connected_weights2 = W_fc2
            b_fc2 = init_bias('biases',[self.num_fc2_neurons])
            h_fc2 = tf.nn.relu(tf.matmul(h_fc1,W_fc2)+b_fc2,name=scope.name)
            activation_summary(h_fc2)

        ## Dropout - Apply dropout before readout layer to reduce overfitting
        with tf.variable_scope('dropout') as scope:
            #keep_prob = tf.placeholder(tf.float32,name="keep_prob")
            h_fc2_drop = tf.nn.dropout(h_fc2, keep_prob,name=scope.name)

        # Readout
        with tf.variable_scope('readout_logits') as scope:
            W_fcr = init_weights('weights',[self.num_fc2_neurons,self.num_classes])
            b_fcr = init_bias('biases',[self.num_classes])
            logits=tf.matmul(h_fc2_drop, W_fcr) + b_fcr
            activation_summary(logits)
        return logits
 




    
    def loss(self,logits, labels):
        #cross_entropy = tf.nn.softmax_cross_entropy_with_logits(logits, labels,
        #                                                        name='cross_entropy')

        #Note: For binary classification I believe
        #sigmoid should work fine in place of softmax in terms of effectiveness
        #Softmax just normalizes over the different classes. If there's only
        #two classes, the values will simply be [-p,p]
        cross_entropy = -tf.reduce_sum(labels*tf.log(tf.clip_by_value(tf.nn.softmax(logits),1e-10,1.0)))
        #To use scalar summary, first argument needs to be a list
        #with same shape as cross_entropy
        #tf.scalar_summary(cross_entropy.op.name, cross_entropy)
        #cross_entropy = -tf.reduce_sum(labels * tf.log(logits), reduction_indices=[1])
        loss = tf.reduce_mean(cross_entropy,
                              name='xentropy_mean')
        activation_summary(loss)
        return loss

    def training(self,loss,learning_rate):
        #Create a scalar summary for loss function
        tf.scalar_summary(loss.op.name, loss)
        # Create the gradient descent optimizer with the given learning rate.
        optimizer = tf.train.GradientDescentOptimizer(learning_rate)
        # Create a variable to track the global step.
        global_step = tf.Variable(0, name='global_step', trainable=False)
        train_op = optimizer.minimize(loss,global_step = global_step)
        return train_op

    def training_adam(self,loss,learning_rate):
        #Create a scalar summary for loss function
        tf.scalar_summary(loss.op.name, loss)
        # Create the gradient descent optimizer with the given learning rate.
        optimizer = tf.train.AdamOptimizer(learning_rate)
        # Create a variable to track the global step.
        global_step = tf.Variable(0, name='global_step', trainable=False)
        train_op = optimizer.minimize(loss,global_step = global_step)
        return train_op


    def logits_to_probs(self,logits):
        return tf.sigmoid(logits)

    def evaluation(self,logits, labels):
        """Evaluate the quality of the logits at predicting the label.
        Args:
        logits: Logits tensor, float - [batch_size, NUM_CLASSES].
        labels: Labels tensor, int32 - [batch_size], with values in the
        range [0, NUM_CLASSES).
        Returns:
        A scalar int32 tensor with the number of examples (out of batch_size)
        that were predicted correctly.
        """
        # For a classifier model, we can use the in_top_k Op.
        # It returns a bool tensor with shape [batch_size] that is true for
        # the examples where the label is in the top k (here k=1)
        # of all logits for that example.
        correct = tf.equal(tf.argmax(logits,1), tf.argmax(labels,1))
        #correct = tf.nn.in_top_k(logits, tf.cast(labels,tf.int32), 1)
        # Return the number of true entries.
        return tf.reduce_sum(tf.cast(correct, tf.int32))





#Convenience functions across classes

def init_weights(name,shape,mean=0,stddev=0.1):
    #Initialize weight variables from a truncated normal distribution
    initial = tf.truncated_normal(shape,mean=0,stddev=stddev)
    return tf.Variable(initial,name=name)

def init_bias(name,shape):
    #Initialize bias variables with a value of 0.1
    initial = tf.constant(0.1, shape=shape)
    return tf.Variable(initial,name=name)


#Initialize weights with decay (only for models with weight decay
def init_weights_decay(name, shape, stddev, wd):
    """Helper to create an initialized Variable with weight decay.
    Note that the Variable is initialized with a truncated normal distribution.
    A weight decay is added only if one is specified.
    Args:
    name: name of the variable
    shape: list of ints
    stddev: standard deviation of a truncated Gaussian
    wd: add L2Loss weight decay multiplied by this float. If None, weight
    decay is not added for this Variable.
    Returns:
    Variable Tensor
    """
    var = _variable_on_cpu(name, shape,
                            tf.truncated_normal_initializer(stddev=stddev))
    if wd is not None:
        weight_decay = tf.mul(tf.nn.l2_loss(var), wd, name='weight_loss')
        tf.add_to_collection('losses', weight_decay)
    return var



def activation_summary(in_op):
    """Helper to create summaries for activations.
    Creates a summary that provides a histogram of activations.
    Creates a summary that measure the sparsity of activations.
    Args:
    x: Tensor
    Returns:
    nothing
    """
    # Remove 'tower_[0-9]/' from the name in case this is a multi-GPU training
    # session. This helps the clarity of presentation on tensorboard.
    #tensor_name = re.sub('%s_[0-9]*/' % TOWER_NAME, '', x.op.name)
    tensor_name = in_op.op.name
    tf.histogram_summary(tensor_name + '/activations', in_op)
    tf.scalar_summary(tensor_name + '/sparsity', tf.nn.zero_fraction(in_op))
