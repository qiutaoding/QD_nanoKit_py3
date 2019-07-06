#!/usr/bin/env python3
#Data 2019/07/06
#Auther: Qiutao DING, PhD candicate
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
    #counter = 0
    trimid_fh = str("trim/trim_%s.list" % barcoder)
    triminfo_fh = str("trimming_info/%s.csv" % barcoder)
    fasta_fh = str("trim/trim_%s.fa" % barcoder)
    output_fh = str("trimmed/trimmed_%s.fa" % barcoder)
    trimminfo = pd.read_csv(triminfo_fh, sep = "\t", index_col = 0)
    record_dict = SeqIO.index(fasta_fh, "fasta")
    with open(trimid_fh, "r") as idlist:
        with open(output_fh,"w") as fw:
            for line in idlist:
                readname = line.rstrip()
                #counter += 1
                #if counter % 500 == 0 : print(str("Processed %i") % (counter))
                readinfo = trimminfo.loc[readname,]
                qstart = readinfo.start
                qstop = readinfo.stop
                qstrand = readinfo.strand
                read = record_dict[readname]
                seq_L = ""
                seq_R = ""
                if qstrand == "+":
                    if len(read) - qstop >= 100 and qstart >= 100:
                        seq_L = str(read.seq)[:qstart]
                        header_L = str(">%s_L len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                        seq_R = str(read.seq)[qstop + 1:]
                        header_R = str(">%s_R len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n")
                    elif len(read) - qstop >= 100: 
                        # right flanking only
                        seq_R = str(read.seq)[qstop + 1:]
                        header_R = str(">%s_R_only len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n")
                    elif qstart >= 100: 
                        # left flanking only
                        seq_L = str(read.seq)[:qstart]
                        header_L = str(">%s_L_only len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                elif qstrand == "-":
                    if len(read) - qstop >= 100 and qstart >= 100: 
                        seq_L = str(read.seq[qstop+1:].reverse_complement())
                        seq_R = str(read.seq[:qstart].reverse_complement())
                        header_L = str(">%s_L len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                        header_R = str(">%s_R len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n")
                    elif len(read) - qstop >= 100: 
                        seq_L = str(read.seq[qstop+1:].reverse_complement())
                        header_L = str(">%s_L_only len:%i strand:%s\n" % (read.id, len(seq_L), qstrand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                    elif qstart >= 100: 
                        seq_R = str(read.seq[:qstart].reverse_complement())
                        header_R = str(">%s_R_only len:%i strand:%s\n" % (read.id, len(seq_R), qstrand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n") 
        fw.close()
    idlist.close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """
    will be added in the furture
    """)
    parser.add_argument("-b", "--barcode", help = "barcode number")
    args = parser.parse_args()
    #function calling
    trimming(fileName(args.barcode))
