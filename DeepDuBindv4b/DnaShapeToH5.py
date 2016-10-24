import numpy as np
import h5py
import os
import glob
import sys


#ext_list = ['.MGW','.Roll','.ProT','.HelT']
ext_list = ['.Roll','.ProT','.HelT']

def main(argv):
    #Note: The fasta file must contain only one entry
    run_dir = "ENCODE_DREAM_LINK/annotations/chromosomes"
    #run_dir = "OLD"
    #run_dir = sys.argv[1]
    #fname = "ENCODE_DREAM_LINK/annotations/chromosomes/chr13.fa.Roll"
    #fname = "fasta_example_one.fa.MGW"
    convert_dnashapes_by_extension(run_dir,'hg19')
    #convert_dnashapes_on_dir(run_dir)

   
def convert_dnashapes_by_extension(run_dir,out_file_base):
    #Create a single HDF5 file for every file extension in ext_list
    os.chdir(run_dir)
    #ext_list = ['MGW','Roll','ProT','HelT']
    
    for ext in ext_list:
        out_file_name = out_file_base+ext+'.h5'
        for fname in glob.glob(("*.fa"+ext)):
            print 'Writing', fname, 'to', out_file_name
            contig_name = os.path.splitext(os.path.splitext(os.path.basename(fname))[0])[0]
            #Contig name example:'chr10'
            with h5py.File(out_file_name,'a') as hf:
                try:
                    hf.create_dataset(contig_name,
                                      data=read_file_to_numpy(fname),
                                      chunks=True)
                except ValueError:
                    print "Likely lack of uniformity in line width in", fname
                    print "Input files must contain only one fasta entry each."
    
            

        
def convert_dnashape(fname):
    ext = os.path.splitext(fname)[1]
    if (ext in ext_list)==False:
        print "Incorrect extension"
        return None  
        
    data_name = os.path.basename(fname)
    print "Converting", data_name, "to h5 format..."
    out_file = fname+".h5"
    with h5py.File(out_file,'w') as hf:
        try:
            hf.create_dataset(data_name,data=read_file_to_numpy(fname),
                                  chunks=True)
        except ValueError:
            print "Likely lack of uniformity in line width in", fname
            print "Input files must contain only one fasta entry each."
           
            
def convert_dnashapes_on_dir(run_dir):
    #Attempts to convert every single file matching the proper extensions
    #to dna shape files
    os.chdir(run_dir)
    #ext_list = ['MGW','Roll','ProT','HelT']

    for ext in ext_list:
        for fname in glob.glob(("*.fa"+ext)):
            convert_dnashape(fname)
                
def read_file_to_numpy(fname):
        np_arr = np.genfromtxt(fname,dtype=None,skip_header=1,skip_footer=1,
                               comments='#',
                               delimiter=',',unpack = True, missing_values=('NA'),
                               filling_values=0.0,autostrip = True)
        last_line = read_last_line_as_numpy(fname)
        #print np.concatenate( (np_arr.transpose().ravel(),last_line),axis=0)
        return np.concatenate( (np_arr.transpose().ravel(),last_line),axis=0)
    
    
    
def read_last_line_as_numpy(fname):
    #Read last line of file as numpy array. Converts 'NA' entries to 0.00
    last_line = read_last_line(fname).split(',')
    for i,item in enumerate(last_line):
        if last_line[i] == 'NA' or last_line[i] == 'NA\n':
            last_line[i] =0.0
    return np.asarray(last_line)
   
def read_last_line(fname):
    #Super fast method to read last line of a file.
    #http://stackoverflow.com/questions/3346430/
    #what-is-the-most-efficient-way-to-get-first-and-last-line-of-a-text-file

    with open(fname,'r') as f:
        f.seek(-2,2) #jump to second to last byte in file
        while f.read(1) != "\n": #until EOL is found
            f.seek (-2,1) #Jump back two bytes (actually one byte)
                          #Note we jump back two because checking
                          #and reading pushes up one
        return f.readline()
                
                
if __name__ == "__main__":
    main(sys.argv)        
