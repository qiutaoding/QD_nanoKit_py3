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
def getIndex(reference):
    """
    Find reference sequence
    make a index
    return mappy alignment result
    """
    if reference:
        reffa = reference
    else:
        reffa = path.join(path.dirname(path.abspath(path.dirname(__file__))),"reference.fa")
    if not path.isfile(reffa):
        logging.error("Could not find reference.fa")
        sys.exit("ERROR: Could not find reference.fa! Programme exit due to reference.fa problem.")
    aligner=mp.Aligner(reffa, preset = "map-ont", best_n = 1, n_threads = 4) #invoke minimap2 API from Li Heng, keep only best alignmen 
    if not aligner:
        logging.error("Failed to load/build index")
        raise Exception("ERROR: failed to load/build index! Programme exit due to mappy problem.")
    return aligner

def align(aligner, reads, output):
    """
    Check reads can be mapped to vector sequence
    """
    counter = 0
    with open(os.path.join(output,"trimmed.fa"),"w") as fw:
        for read in SeqIO.parse(reads, "fasta"):  
            try:
                alignInfo = next(aligner.map(str(read.seq)))
                counter +=1
                seq_L=''
                seq_R=''
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
            except StopIteration:
                print(read.format("fasta"), end='')
    fw.close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Remove reference sequence in read and save the trimmed to output directory, and filter sequence without reference sequence to standard output",\
    epilog="Quick start: FosmidSequencingRemoveVector.py -r reference.fa -i input.fa -o output_directory > sequence_without_vector.fa")
    parser.add_argument("-r", "--reference", help = "reference fasta")
    parser.add_argument("-i", "--input", help = "input fasta for trimming")
    parser.add_argument("-o", "--output", help = "output directory to save trimming information")
    args = parser.parse_args()
    #function calling
    aligner = getIndex(args.reference)
    align(aligner, args.input, args.output)
    print("Job finished.")