####### Script to read strain and sample names for a file and copy appropriate reads from the directory in output directory ###
import sys
import read_subdirectories_names
import glob

strain_name_file=sys.argv[1]    #### tab-delimited strain and sample name of read files ###
read_input_dir=sys.argv[2]       ####  directory path to read from ######
read_output_dir=sys.argv[3]     #### directory pathy to write into #####

#### read input directories ####

list_subdirectories=read_subdirectories_names.get_immediate_subdirectories(read_input_dir)

##### open and read strain and sample file ####

with open(strain_name_file,'r') as list_strain_name:

    for line in list_strain_name:
    
        #print line
        column_list=line.split('\t')
        
        print column_list
    
        strain_name=column_list[0]
        dir_name=column_list[1]
        
        print strain_name
        print dir_name
        
        if glob.glob(read_input_dir+'*/dir_name*'):
            
            print 1
            
        else:
            
            print 0
        
     
            
            
       
            
        
        


