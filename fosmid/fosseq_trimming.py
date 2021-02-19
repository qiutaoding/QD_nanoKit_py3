#!/usr/bin/env python3
#Data 2019/07/14
#Auther: Qiutao DING, PhD candidate
#Department of Biology, Hong Kong Baptist University, Hong Kong, China
#Contact: qiutaod@gmail.com 
import os
import pandas as pd
from Bio import SeqIO
import argparse

def fileName(j=3):
    if int(j) <10:
        i = "barcode0"+str(j)
    else:
        i = "barcode"+str(j)
    return(i)

def trimming(barcoder):
    triminfo_fh = str("blastn_backbone_info_filtered/%s.csv" % barcoder)
    fasta_fh = str("read_with_backbone/%s.fa" % barcoder)
    output_fh = str("read_trimmed/trimmed_%s.fa" % barcoder)
    trimminfo = pd.read_csv(triminfo_fh, sep = "\t", index_col = 0,header=None,names=["start","end","strand","alignlength"])
    record_dict = SeqIO.index(fasta_fh, "fasta")
    with open(output_fh,"w") as fw:
        for readid in trimminfo.index:
            readinfo = trimminfo.loc[readid,]
            qstart = readinfo.start
            qstop = readinfo.end
            qstrand = readinfo.strand
            read = record_dict[readid]
            seq_L = ""
            seq_R = ""
            if qstrand == "+":
                if len(read) - qstop >= 100 and qstart >= 100:
                    seq_L = str(read.seq)[:qstart-1]
                    header_L = str(">%s_L len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                    fw.write(header_L)
                    fw.write(seq_L+"\n")
                    seq_R = str(read.seq)[qstop :]
                    header_R = str(">%s_R len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                    fw.write(header_R)
                    fw.write(seq_R+"\n")
                elif len(read) - qstop >= 100: 
                    # right flanking only
                    seq_R = str(read.seq)[qstop :]
                    header_R = str(">%s_R len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                    fw.write(header_R)
                    fw.write(seq_R+"\n")
                elif qstart >= 100: 
                    # left flanking only
                    seq_L = str(read.seq)[:qstart-1]
                    header_L = str(">%s_L len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                    fw.write(header_L)
                    fw.write(seq_L+"\n")
            elif qstrand == "-":
                if len(read) - qstop >= 100 and qstart >= 100: 
                    seq_L = str(read.seq[qstop:].reverse_complement())
                    seq_R = str(read.seq[:qstart-1].reverse_complement())
                    header_L = str(">%s_L len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                    fw.write(header_L)
                    fw.write(seq_L+"\n")
                    header_R = str(">%s_R len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                    fw.write(header_R)
                    fw.write(seq_R+"\n")
                elif len(read) - qstop >= 100: 
                    seq_L = str(read.seq[qstop:].reverse_complement())
                    header_L = str(">%s_L len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                    fw.write(header_L)
                    fw.write(seq_L+"\n")
                elif qstart >= 100: 
                    seq_R = str(read.seq[:qstart-1].reverse_complement())
                    header_R = str(">%s_R len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                    fw.write(header_R)
                    fw.write(seq_R+"\n") 
    fw.close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """
    will be added in the furture
    """)
    parser.add_argument("-b", "--barcode", help = "barcode number")
    args = parser.parse_args()
    #function calling
    trimming(fileName(args.barcode))
    print("Finished",fileName(args.barcode))
