#!/usr/bin/env python3
#Data 2019/06/24
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

def align(aligner, reads, output):
    """
    Check reads can be mapped to vector sequence
    """
    counter = 0
    with open(os.path.join(output),"w") as fw:
        for read in SeqIO.parse(reads, "fasta"):  
            try:
                alignInfo = next(aligner.map(str(read.seq)))
                counter +=1
                seq_L=''
                seq_R=''
                if alignInfo.strand == 1:
                    if len(read) - alignInfo.q_en >= 100 and alignInfo.q_st >= 100: # if positive strand mapping to pSMART backbone with both genomic flankings longer than 100 bp
                        seq_L = str(read.seq)[:int(alignInfo.q_st)-1]
                        seq_R = str(read.seq)[int(alignInfo.q_en):]
                        header_L = str(">%s_L len:%i strand:%i\n" % (read.id, len(seq_L), alignInfo.strand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                        header_R = str(">%s_R len:%i strand:%i\n" % (read.id, len(seq_R), alignInfo.strand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n")
                    elif len(read) - alignInfo.q_en >= 100: # if right flanking is longer than 100 bp ,with shorted than 100 bp left flanking or no left flanking sequence
                        seq_R = str(read.seq)[int(alignInfo.q_en):]
                        header_R = str(">%s_R_only len:%i strand:%i\n" % (read.id, len(seq_R), alignInfo.strand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n")
                    elif alignInfo.q_st >= 100: # if left flanking is longer than 100 bp ,with shorted than 100 bp right flanking or no left flanking sequence
                        seq_L = str(read.seq)[:int(alignInfo.q_st)-1]
                        header_L = str(">%s_L_only len:%i strand:%i\n" % (read.id, len(seq_L), alignInfo.strand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                elif alignInfo.strand == -1: #if read is reverse complementary to vector sequence
                    if len(read) - alignInfo.q_en >= 100 and alignInfo.q_st >= 100: # if read mapping to pSMART backbonenegative strand with both genomic flankings longer than 100 bp
                        seq_L = str(read.seq.reverse_complement())[:len(read)-int(alignInfo.q_en)]
                        seq_R = str(read.seq.reverse_complement())[len(read)-int(alignInfo.q_st)+1:]
                        header_L = str(">%s_L len:%i strand:%i\n" % (read.id, len(seq_L), alignInfo.strand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                        header_R = str(">%s_R len:%i strand:%i\n" % (read.id, len(seq_R), alignInfo.strand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n")
                    elif len(read) - alignInfo.q_en >= 100: # if left flanking is longer than 100 bp ,with shorted than 100 bp left flanking or no left flanking sequence
                        seq_L = str(read.seq.reverse_complement())[:len(read)-int(alignInfo.q_en)]
                        header_L = str(">%s_L len:%i strand:%i\n" % (read.id, len(seq_L), alignInfo.strand))
                        fw.write(header_L)
                        fw.write(seq_L+"\n")
                    elif alignInfo.q_st >= 100: # if left flanking is longer than 100 bp ,with shorted than 100 bp right flanking or no left flanking sequence
                        seq_R = str(read.seq.reverse_complement())[len(read)-int(alignInfo.q_st)+1:]
                        header_R = str(">%s_R len:%i strand:%i\n" % (read.id, len(seq_R), alignInfo.strand))
                        fw.write(header_R)
                        fw.write(seq_R+"\n")      
            except StopIteration:
                print(read.format("fasta"), end='')
    fw.close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Remove reference sequence in read and save the trimmed to output directory, and filter sequence without reference sequence to standard output",\
    epilog="Quick start: FosmidSequencingRemoveVector.py -r fosmid_sample/pSMART_BAC.fa -i fosmid_sample/sample03.fa -o fosmid_sample/result/trimming.fa > fosmid_sample/result/non_vector_reads.fa")
    parser.add_argument("-r", "--reference", help = "reference fasta")
    parser.add_argument("-i", "--input", help = "input fasta for trimming")
    parser.add_argument("-o", "--output", help = "output file to save trimming information")
    parser.add_argument("-t", "--thread", help = "threads N want to use during alignment")
    args = parser.parse_args()
    #function calling
    aligner = getIndex(args.reference, args.thread)
    align(aligner, args.input, args.output)