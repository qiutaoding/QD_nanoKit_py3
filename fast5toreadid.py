#!/usr/bin/env python3
#Data 2019/06/20
#Auther: Qiutao DING, PhD candicate
#Department of Biology, Hong Kong Baptist University, Hong Kong, China
#Contact: qiutaod@gmail.com 
import os
import h5py
import argparse
import csv
# List all files in a directory using os.listdir
# This function convert nanopore fast5 strand name to basecalled fastq read name
def fast5scan(basepath,savepath):
#scan directory in input directory
    with os.scandir(basepath) as subdir:
        fold_counter=0
        dict_inti={}
    #reading subdirectory
        for entry in subdir:
            #reading files under subdirectory
            temp_dir=os.path.join(entry)
            with os.scandir(temp_dir) as filelist:
                for subfile in filelist:
                    if subfile.is_file():
                        FileHandle=h5py.File(subfile,"r")
                        read_name=""
                        read_group=""
                        #reading fast5 subgroup
                        for readgroup in FileHandle["Raw/Reads/"]:
                            read_group=readgroup
                            #extract read id from attributes
                        for CharASCII in FileHandle[str("Raw/Reads/"+read_group)].attrs["read_id"]:
                            read_name+=chr(CharASCII)
                        #add key and value to folder dict
                        dict_inti[str(subfile)[11:-8]]=str(read_name)
            if fold_counter % 1 ==0: #print progress
                print("Processing folder",fold_counter)
            with open(os.path.join(savepath,"fast5toreadname_"+str(fold_counter)+".csv"),"w") as fw:
                file_writer=csv.writer(fw)
                filednames=['fast5','fastq_read_name']
                file_writer.writerow(filednames)
                for key in dict_inti.keys():
                    file_writer.writerow([key,dict_inti[key]])
            fw.close()
            dict_inti={}
            fold_counter+=1
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """
    Quick start: 
    read_fast5_folder_build_read_id_table.py -i fast5_directory/ -s saving_directory/
    """)
    parser.add_argument("-i", "--input", help = "input directory")
    parser.add_argument("-s", "--save", help = "saving directory")
    args = parser.parse_args()
    #function calling
    fast5scan(args.input,args.save)