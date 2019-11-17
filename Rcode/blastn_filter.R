setwd("F:/nanopore_data/JUfosmid/seperation")
library(data.table)
library(readr)
library(dplyr)
require(intervals)
header = c('qname', 'rname','pident','length','mismatch','gapopen','qstart','qend','rstart','rend','bit','score')
j=3

fileName = paste0('map_info/barcode',formatC(j, width = 2, flag = '0'),'.csv')
saving = paste0('map_info/barcode',formatC(j, width = 2, flag = '0'),'_info.csv')
df <- read.csv(file = fileName, header = F, sep = '\t')
colnames(df)<-header
df_ex <-subset(df, select = c(1,2,4,7,8,9,10))
rm(df)

readlist_fileName = paste0('map_info/barcode',formatC(j, width = 2, flag = '0'),'.list')
readInList = file(readlist_fileName, "r")
readList <-readInList
read_num= 0
blastn_df <-df_ex[FALSE,]
while (TRUE) {
  line = readLines(readList, n=1)
  if(length(line)==0){break}else{
    read_name = paste0(line)
    read_num = read_num + 1
    if (read_num %% 100 == 0){print(paste("Working barcode",formatC(j, width = 2, flag = '0')," read:", read_num))}
    df_tmp <- subset(df_ex, qname == read_name)
    filter_df <-df_tmp[FALSE,]
    filter_df$qname = read_name
  }
}
close(readInList)




read_name = "1d2df738-538f-47a2-bf09-f6b6e4e5a4d6"
df_tmp <- subset(df_ex, qname == read_name)

#df_tmp <-df_tmp[order(df_tmp$length),]
#df_filter <- rbind(df_filter,df_tmp[nrow(df_tmp),])
View(df_filter)


blastn_df = data.frame(readname = character(),
                       start = integer(), 
                       end = integer(), 
                       strand = character(),
                       stringsAsFactors=FALSE)
colnames(blastn_df) = c("readname", "start", "end", "strand")


read_map_union <- as.data.frame(interval_union(Intervals(df_tmp[,4:5])))
read_map_union_df <-cbind(qname = read_name, read_map_union, strand = NA)
colnames(read_map_union_df)<-c("qname","qstart","qend","strand")
for (n in seq(1:nrow(read_map_union_df))){
  tmp_strand_df <- subset(df_tmp, qstart == read_map_union_df[n,]$qstart )
  if (tmp_strand_df[1,]$rend - tmp_strand_df[1,]$rstart > 0){
    read_map_union_df[n,]$strand <- "+"
  }else{
    read_map_union_df[n,]$strand <- "-"
  }
  
  blastn_df <-rbind(blastn_df, read_map_union_df)
}









