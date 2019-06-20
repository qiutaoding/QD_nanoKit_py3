#!/usr/bin/env python3
#Data 2019/06/18
#Auther: Qiutao DING, PhD candicate
#Department of Biology, Hong Kong Baptist University, Hong Kong, China
#Contact: qiutaod@gmail.com 

from Bio import SeqIO
import argparse
import os
import sys
def FastaLength(FileName):
    """
    Target a NNNN.fasta or NNNN.fa file
    This script is designed for small genome, i.e. Nematode (C. elegans)
    Return list of sequence length
    """
    sequence_length_list = []
    flush_n=1
    print("Reading fasta file")
    FileHandle = open(FileName, "r")
    for sequence in SeqIO.parse(FileHandle, "fasta"):
        sys.stdout.write(str(flush_n) + "-")
        sys.stdout.flush()        
        sequence_length_list.append(len(sequence))
        flush_n += 1
    FileHandle.close
    sys.stdout.write(str(flush_n))
    print("")
    return sequence_length_list

def N50(length_list):
    """
    get N50 read length list from the fasta file
    """
    #calculate N50 value
    N50_boundary = sum(length_list)*0.5
    #sort read length by length value
    length_list.sort()
    counter = 0
    for length_info in length_list:
        counter += length_info
        #stop adding length information after over N50 value
        if counter >=N50_boundary:
            break
    print("")
    return length_info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Quick start: N50_calculation.py NNNN.fa")
    parser.add_argument("-i", "--input", help = "input file")
    args = parser.parse_args()
    #function calling
    print("N50 is {:,} bp".format(N50(FastaLength(args.input))))