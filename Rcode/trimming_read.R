setwd("F:/nanopore_data/JUfosmid/seperation")
library(data.table)
library(readr)
library(dplyr)
header = c('qname', 'qlen','qstart','qend','strand','rname','rlen','rstart','rend')

j=3
###################filter df for further read process############################
for (j in seq(1,57)){
  Align_df_fileName = paste0('align_bac/bac_barcode',formatC(j, width = 2, flag = '0'),'.paf')
  Align_df <- read.csv(Align_df_fileName,sep = "\t", header = F)
  Align_df_trim <- subset(Align_df, select = c(seq(1,9)))
  rm(Align_df)
  colnames(Align_df_trim) <-header
  #View(Align_df_trim)
  readlist_fileName = paste0('trim/trim_barcode',formatC(j, width = 2, flag = '0'),'.list')
  readInList = file(readlist_fileName, "r")
  readList <-readInList
  read_num = 0
  fos_trimming_df = data.frame("readName" = character(),
                               "length" = integer(),
                               "start" = integer(),
                               "stop" = integer(), 
                               "strand" = character(),
                               stringsAsFactors=FALSE)
  while (TRUE) {
    line = readLines(readList, n=1)
    if(length(line)==0){break} 
    read_name = paste0(line)
    read_num = read_num + 1
    working_df = subset(Align_df_trim, qname == read_name)
    if (table(working_df$strand)["+"] > table(working_df$strand)["-"]){
      add_info = c(read_name ,
                   max(working_df$qlen),
                   min(working_df$qstart),
                   max(working_df$qend),
                   "+" )
    }else{
      add_info = c(read_name , 
                   max(working_df$qlen),
                   min(working_df$qstart), 
                   max(working_df$qend), 
                   "-" )
    }
    fos_trimming_df[nrow(fos_trimming_df)+1 , ]= add_info
  }
  close(readInList)
  
  print(paste("barcode",j,"read count:", read_num))
  saving = paste0('trimming_info/barcode',formatC(j, width = 2, flag = '0'),'.csv')
  write.table(fos_trimming_df, file = saving, sep = '\t',row.names=FALSE, col.names = TRUE)
}




