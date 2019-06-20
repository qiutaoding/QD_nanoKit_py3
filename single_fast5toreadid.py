#!/usr/bin/env python3
#Data 2019/06/20
#Auther: Qiutao DING, PhD candicate
#Department of Biology, Hong Kong Baptist University, Hong Kong, China
#Contact: qiutaod@gmail.com 
#!/usr/bin/env python3
import os
import h5py
import argparse
# List all files in a directory using os.listdir
# This function convert nanopore fast5 strand name to basecalled fastq read name
def fast5scan(filepath):
    FileHandle=h5py.File(filepath,"r")
    read_name=""
    #reading fast5 subgroup
    for readgroup in FileHandle["Raw/Reads/"]:
        read_group=readgroup
        #extract read id from attributes
    for CharASCII in FileHandle[str("Raw/Reads/"+read_group)].attrs["read_id"]:
        read_name+=chr(CharASCII)
    return(read_name)

            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """
    Quick start: 
    single_fast5toreadid.py -i fast5_test/0/GXB01298_20190601_FAK27254_GA10000_sequencing_run_juFosSP1_24_85131_read_4_ch_7_strand.fast5
    """,
    epilog="Result is shown in terminal")
    parser.add_argument("-i", "--input", help = "input testing fast5 file")
    args = parser.parse_args()
    #function calling
    print(fast5scan(args.input))