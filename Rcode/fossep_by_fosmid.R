setwd("F:/nanopore_data/JUfosmid/seperation")
library(data.table)
library(readr)
library(dplyr)

#Primary filtering trimmed read overlapping information from minimap2 
#by 1. removing self-overlap ,and 2. reverse strand overlapping.
header = c('qname', 'qlen','qstart','qend','strand','rname','rlen','rstart','rend')
for (i in seq(1,57)){
  fileName = paste0('flanking_align/trim_barcode',formatC(i, width = 2, flag = '0'),'.paf')
  saving = paste0('trim_paf2csv/filter_barcode',formatC(i, width = 2, flag = '0'),'.csv')
  df <- read.csv(file = fileName, header = F, sep = '\t')
  df_ex <-subset(df, select = c(seq(1,9)))
  rm(df)
  #in some cases, left flanking or right flanking is long enough to align to other flanking sequence
  df_ex2 <- subset(df_ex, df_ex$V1 != paste0(df_ex$V6) & df_ex$V5 == "+")
  rm(df_ex)
  colnames(df_ex2) <- header
  df_ex2$qper <- with(df_ex2, abs(qend-qstart)/qlen )
  df_ex2$rper <- with(df_ex2, abs(rend-rstart)/rlen )
  write.table(df_ex2, file = saving, sep = '\t',row.names=FALSE)
  print(paste("Complete", formatC(i , width = 2 , flag = '0')))
}


#stopped here at 2019/07/06 17:12
########################main content###############################
##build functions
#building fucntions to extract read name
position_set = 400
cutoff_set = 0.70
extractAlign <- function(read_name){
  single_foslist_initial <-list()
  df_read_initial<-rbind(df_filter[df_filter$qname %like% read_name, ],
                         df_filter[df_filter$rname %like% read_name, ])
  tmp1 <- df_read_initial %>% filter(substr(qname , start = 38 , stop =38)=="R",
                                     substr(rname , start = 38 , stop =38)=="R", 
                                     qstart < position_set , rstart < position_set,
                                     qper > cutoff_set | rper > cutoff_set)
  tmp2 <- df_read_initial %>% filter(substr(qname , start = 38 , stop =38)=="L",
                                     substr(rname , start = 38 , stop =38)=="L", 
                                     qlen-qend < position_set , rlen-rend < position_set,
                                     qper > cutoff_set | rper > cutoff_set)
  
  tmp = rbind(tmp1,tmp2)
  if (nrow(tmp) == 0) {return(single_foslist_initial)}else{
    #make mapped fosmid flankings read name as a unique list, extract read name excluding _L or _R
    single_foslist_initial <-unique(c(substr(tmp$qname,start = 1 , stop = 36),
                                      substr(tmp$rname,start = 1 , stop = 36)))
    return(single_foslist_initial)}}

looping_foslist <-function(initial_foslist){
  final_foslist = list()
  #pass input list to return list to add new element
  if (length(initial_foslist)==0){
    return(initial_foslist)
  }else{
    final_foslist = c(final_foslist,read_name)
    for (read_screening in initial_foslist){
      if (read_screening == read_name){
        next
      }else{
        single_foslist_screen <- extractAlign(read_screening)
        for (screen_element in single_foslist_screen){
          if (screen_element %in% unlist(final_foslist) |
              screen_element %in% fos_sep_df$read) {next
          }else{
            final_foslist <-c(final_foslist, screen_element)}
        }}}}
  return(final_foslist)
}

#main()
for (j in seq(1,57)){
  FilereadInAlign_fileName = paste0('trim_paf2csv/filter_barcode',formatC(j, width = 2, flag = '0'),'.csv')
  FilereadInAlign <- read.csv(FilereadInAlign_fileName,sep = "\t", header = T)
  df_filter <- subset(FilereadInAlign, strand == '+' )
  #save ram
  rm(FilereadInAlign)
  readlist_fileName = paste0('trimmed/flanking_barcode',formatC(j, width = 2, flag = '0'),'.list')
  readInList = file(readlist_fileName, "r")
  readList <-readInList
  fosmid_counter = 0
  fosmid_counter_single = 0
  #initialize a empty list for fosmid seperation
  fos_sep_df = data.frame(fos_number = integer(), read=character(), stringsAsFactors=FALSE)
  while (TRUE) {
    line = readLines(readList, n=1)
    if(length(line)==0){break}else{
      read_name = paste0(line)
      if (read_name %in% fos_sep_df$read){next}else{
        #print progress
        uniq_list <- looping_foslist(extractAlign(read_name))
        if (length(uniq_list) == 0 ){
          #add read without best alignment to list
          fosmid_counter_single = fosmid_counter_single + 1
          fosmid_counter = fosmid_counter + 1
          fos_sep_df[nrow(fos_sep_df)+1, ] =c(fosmid_counter, read_name)
        }else{
          fosmid_counter = fosmid_counter + 1
          for (element in uniq_list){fos_sep_df[nrow(fos_sep_df)+1, ] =c(fosmid_counter, element)}}}}}
  close(readInList)
  print(paste0("Finished barcode",formatC(j, width = 2, flag = '0')," with ",
               fosmid_counter, " fosmid with ",fosmid_counter_single, " singleton."  ))
  #save list to csv
  write.table(data.frame(fos_sep_df),
              paste0('fos_sep/sample',formatC(j,width=2,flag='0'),'.csv'),
              col.names=T,sep='\t')
}





