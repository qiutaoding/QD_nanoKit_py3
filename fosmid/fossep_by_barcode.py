#!/usr/bin/env python3
#Data 2019/07/06
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

#extract df which qname or rname containing read name
def extractAlign(readname):
    single_foslist_initial = list()
    read_align_q = pd.DataFrame()
    read_align_r = pd.DataFrame()
    read_align_q = flanking_align[flanking_align["qname"].str.contains(readname)]
    read_align_r = flanking_align[flanking_align["rname"].str.contains(readname)]
    df_read_initial = read_align_q.append(read_align_r)
    #filter fosmid flanking alignment are parallelly aligned
    df_tmp = df_read_initial[((df_read_initial["qname" and "rname"].str.contains("R")) & 
                    (df_read_initial["qstart"] < position_set) & 
                    (df_read_initial["rstart"] < position_set) & 
                    ((df_read_initial["qper"] > cutoff_set)|(df_read_initial["rper"] > cutoff_set))) |
                    ((df_read_initial["qname" and "rname"].str.contains("L")) &
                     (df_read_initial["qlen"] - df_read_initial["qend"] < position_set) &
                     (df_read_initial["rlen"] - df_read_initial["rend"]< position_set) &
                     ((df_read_initial["qper"] > cutoff_set)|(df_read_initial["rper"] > cutoff_set)))]
    #get the initial fosmid sepeartion read list
    if len(df_tmp.index) == 0:
        return(single_foslist_initial)
    else:
        single_fos_list = list(set(list(df_tmp["qname"].str.slice(0,36))+list(df_tmp["rname"].str.slice(0,36))))
        return(single_fos_list)
#looping function for read name in the initial list
def looping_foslist(single_fos_list):
    final_foslist = list()
    if len(single_fos_list) != 0:
        final_foslist = single_fos_list
        for read_screening in single_fos_list:
            if read_screening == readname:
                continue
            else:
                single_foslist_screen = extractAlign(read_screening)
                if len(single_foslist_screen) == 0:
                    continue
                else:
                    for screen_element in single_foslist_screen:
                        if screen_element in final_foslist :
                            continue
                        elif screen_element in fos_df[["readname"]]:
                            continue
                        elif len(screen_element)==0:
                            continue
                        else:
                            final_foslist.append(screen_element)
        return(final_foslist)
    else:
        return(single_fos_list)
#input barcode write fosmid table
def barcodeProcess(barcoder):
    #initialize a counter
    counter = 0
    read_counter = 0
    #file name
    readlist_fh = str("trimmed/flanking_%s.list" % barcoder)
    output_fh = str("python_fos_sep/fos_%s.csv" % barcoder)
    flanking_align_fh = str("trim_paf2csv/filter_%s.csv" % barcoder)
    #input df
    global flanking_align
    flanking_align = pd.read_csv(flanking_align_fh, sep = "\t")
    #make an empty df to add data
    global fos_df
    fos_df = pd.DataFrame(columns=["fosmid", "readname"])
    with open(readlist_fh, "r") as idlist: 
        for line in idlist:
            global readname
            readname = line.rstrip()
            read_counter += 1
            if read_counter % 100 == 0:
                print(str("Complete %s reads."% read_counter))
            if readname in fos_df[["readname"]]:
                continue
            else:
                uniq_list = list()
                uniq_list = looping_foslist(extractAlign(readname))
                if len(uniq_list) == 0:
                    if readname in fos_df[["readname"]]:
                        continue
                    else:
                        counter +=1
                        fos_df = fos_df.append(pd.DataFrame([[counter,readname]],columns=["fosmid", "readname"]))
                else:
                    counter +=1
                    for list_element in uniq_list:
                        fos_df = fos_df.append(pd.DataFrame([[counter,list_element]],columns=["fosmid", "readname"]))
        fos_df.to_csv(output_fh, sep = '\t', index=False)
    idlist.close()
    print(barcoder, "seperate", counter, "fosmid")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = """
    will be added in the furture
    """)
    parser.add_argument("-b", "--barcode", help = "barcode number")
    args = parser.parse_args()
    #function calling
    #set cutoff parameter
    position_set = 400
    cutoff_set = 0.70
    barcoder=fileName(args.barcode)
    barcodeProcess(barcoder)
