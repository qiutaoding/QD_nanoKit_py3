# -*- coding: utf-8 -*-

#define FASTA file 
#1st argment is designed for extracted read name. 2nd is for the sequence in FASTA file
#this readFa function is designed for single fasta sequence parsing
def readFa(readName,readFile):
    filename=str(readFile)+".fa"
    headerInfo="" #sequence header information to return
    sequence=""  #DNA sequences to return
    try:
        with open(filename,"r") as fh:
            for line in fh:
                if line.startswith(">"+readName):
                    headerInfo=line.rstrip()
                    for line in fh:
                        if line.startswith(">"): #if read name found in FASTA, break reading file
                            break
                        sequence+=line.strip()
    except:
        print("ERROR:",filename,"is not found") #Error message if FASTA  is not found 
    else:
        if len(sequence)==0:
            print(readName,"is not found") #Message, if read name is not in the input FASTA
    return(headerInfo,sequence) #return read sequence