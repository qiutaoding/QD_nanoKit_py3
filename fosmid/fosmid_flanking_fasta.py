#!/usr/bin/env python3
#Data 2019/07/13
#Auther: Qiutao DING, PhD candidate
#Department of Biology, Hong Kong Baptist University, Hong Kong, China
#Contact: qiutaod@gmail.com 
import argparse
import pandas as pd
from Bio import SeqIO
def fileName(j=3):
    if j <10:
        i = "barcode0"+str(j)
    else:
        i = "barcode"+str(j)
    return(i)
def sampleName(j=3):
    if j <10:
        i = "sample0"+str(j)
    else:
        i = "sample"+str(j)
    return(i)
def check_read_R(readname):
    if str(readname+"_R_only") in record_dict:
        return (str(readname+"_R_only"))
    elif str(readname+"_R") in record_dict:
        return (str(readname+"_R"))
def check_read_L(readname):    
    if str(readname+"_L_only") in record_dict:
        return (str(readname+"_L_only"))
    elif str(readname+"_L") in record_dict:
        return (str(readname+"_L"))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """
    will be added in the furture
    """)
    parser.add_argument("-b", "--barcode", help = "barcode number")
    args = parser.parse_args()
    barcode = int(args.barcode)
    input_fh = str("trimmed/trimmed_%s.fa" %fileName(barcode))
    flanking_info_fh = str("fosmid_flanking_infocsv/%s.csv" %sampleName(barcode))
    output_L_fh= str("fosmid_flanking_fa/L_%s.fa" %fileName(barcode))
    output_R_fh= str("fosmid_flanking_fa/R_%s.fa" %fileName(barcode))
    global record_dict
    record_dict = SeqIO.index(input_fh, "fasta")
    flanking_info = pd.read_csv(flanking_info_fh, sep = "\t", index_col = 0)
    with open(output_L_fh,"w") as fw_L,open(output_R_fh,"w") as fw_R :
        for index in flanking_info.index:
            if flanking_info["left_flanking"][index] != "empty":            
                header_L = str(">%s_%s_L\n" %(fileName(barcode),index))
                seq_L = str(record_dict[check_read_L(flanking_info["left_flanking"][index])].seq)
                fw_L.write(header_L)
                fw_L.write(seq_L+"\n")
            if flanking_info["right_flanking"][index] != "empty":            
                header_R = str(">%s_%s_R\n" %(fileName(barcode),index))
                seq_R = str(record_dict[check_read_R(flanking_info["right_flanking"][index])].seq)
                fw_R.write(header_R)
                fw_R.write(seq_R+"\n")            
    fw_L.close()
    fw_R.close()
    print(str("=====%s finished=====\n" % (fileName(barcode))))