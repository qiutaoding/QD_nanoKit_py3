#!/usr/bin/env python3
#Data 2019/06/25
#Auther: Qiutao DING, PhD candicate
#Department of Biology, Hong Kong Baptist University, Hong Kong, China
#Contact: qiutaod@gmail.com 
#		  16483960@life.hkbu.edu.hk

import argparse
import os
from os import path
import sys
import mappy as mp
from Bio import SeqIO
def getIndex(reference, thread):
    """
    Find reference sequence
    make a index
    return mappy alignment result
    default only keep one best alignment
    default using 2 threads
    """
    if reference:
        reffa = reference
    else:
        reffa = path.join(path.dirname(path.abspath(path.dirname(__file__))),"reference.fa")
    if not path.isfile(reffa):
        logging.error("Could not find reference.fa")
        sys.exit("ERROR: Could not find reference.fa! Programme exit due to reference.fa problem.")
    if thread is None: #check wether thread arugment is specified, default using 2 threads
        thread = 2
    aligner=mp.Aligner(reffa, preset = "map-ont", best_n = 1, n_threads = int(thread)) #invoke minimap2 API from Li Heng, keep only best alignmen 
    if not aligner:
        logging.error("Failed to load/build index")
        raise Exception("ERROR: failed to load/build index! Programme exit due to mappy problem.")
    return aligner

def align(aligner, reads):
    """
    Check reads can be mapped to vector sequence
    """
    counter = 0
    for read in SeqIO.parse(reads, "fasta"):  
        try:
            alignInfo = next(aligner.map(str(read.seq)))
            print(alignInfo)  
        except StopIteration:
            print(read.format("fasta"), end='')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "minimap2 API Mappy function test, print mapping result",\
    epilog="Quick start: mappy_test.py -r fosmid_sample/pSMART_BAC.fa -i mappp_test_sample/44b73770-d1c6-40f4-979b-8e069e21c2a8.fa")
    parser.add_argument("-r", "--reference", help = "reference fasta")
    parser.add_argument("-i", "--input", help = "input fasta for trimming")
    parser.add_argument("-t", "--thread", help = "threads N want to use during alignment")
    args = parser.parse_args()
    #function calling
    aligner = getIndex(args.reference, args.thread)
    align(aligner, args.input)