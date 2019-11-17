setwd("F:/nanopore_data/JUfosmid/seperation")
#rewrite paf from minimap2 generated alignment information
fileseq <- seq(1,48)
header = c('qname', 'qlen','qstart','qend','strand','rname','rlen','rstart','rend','match','alignlength')
for (i in fileseq){
  fileName = paste0('flanking_align/trim_barcode',formatC(i, width = 2, flag = '0'),'.paf')
  saving = paste0('paf2csv/filter_barcode',formatC(i, width = 2, flag = '0'),'.csv')
  df <- read.csv(file = fileName, header = F, sep = '\t')
  df_ex <-subset(df, select = c(seq(1,11)))
  #in some cases, left flanking or right flanking is long enough to align to other flanking sequence
  df_ex2 <- subset(df_ex, df_ex$V1 != paste0(df_ex$V6))
  colnames(df_ex2) <-header
  df_ex2$qper <- with(df_ex2, abs(qend-qstart)/qlen )
  df_ex2$rper <- with(df_ex2, abs(rend-rstart)/rlen )
  write.table(df_ex2, file = saving, sep = '\t',row.names=FALSE)
  print(paste("Complete", formatC(i , width = 2 , flag = '0')))
}

#convert alignment information into fosmid physical location
#read barcoded reads with backbone sequence, extract their read name
library(data.table)
library(readr)

#fosmid list with at least cutoff number of reads
fos_sep_df = data.frame(barcode = character(), list = integer(), cutoff = integer(),stringsAsFactors=FALSE)
colnames(fos_sep_df) = c("barcode", "list", "list_cutoff")
cutoff = 5
for (j in fileseq){
  #read in read_id file name, for input read name one by one
  readlist_fileName = paste0('trim/trim_barcode',formatC(j, width = 2, flag = '0'),'.list')
  #read in rewrite alignment file for subset and data filtering
  readInAlign_fileName = paste0('paf2csv/filter_barcode',formatC(j, width = 2, flag = '0'),'.csv')
  #read in read_id
  readInList = file(readlist_fileName, "r")
  #read in alignment table
  readInAlign <- read_delim(readInAlign_fileName,"\t", escape_double = FALSE, trim_ws = TRUE)
  #pass read in read_id to readList
  readList <-readInList
  #read number counter
  read_num = 0
  #initilize a counter for fosmid seperation
  fosmid_counter = 0
  #initialize a empty list for fosmid seperation
  fos_list = list()
  while (TRUE) {
    line = readLines(readList, n=1) #read line by line
    if(length(line)==0){break} #stop loop if to the end of file
    read_name = paste0(line) #paste line as read name
    #check if read_name in exsiting fosmid collection
    if (read_name %in% fos_list){
      #need polish condition and process
      next} #go to next read name, if found read name exists in the fosmid collection
    read_num = read_num + 1
    #print progress
    if ( read_num%%50 == 0){print(paste("Working barocode", j ,"read:", read_num))}
    readInFilter <- subset(readInAlign, strand == '+' & (qper >0.80 | rper >0.80))
    tmp1 <- readInFilter[readInFilter$qname %like% read_name, ]
    tmp2 <- readInFilter[readInFilter$rname %like% read_name, ]
    #move to next read id if no alignment found in the file
    if(nrow(tmp1) ==0 & nrow(tmp2) ==0){next}
    #rbind two subset containing read id 
    tmp_merge <- rbind(tmp1, tmp2) 
    #initialize vector
    uniq_list_1 = c()
    uniq_list_2 = c()
    #vector construct
    if(nrow(tmp1) >= 1){
      for (readName in CJ(tmp_merge$qname,  unique = T)){uniq_list_1 <- readName}
      }
    if(nrow(tmp2) >= 1){
      for (readName in CJ(tmp_merge$rname,  unique = T)){uniq_list_2 <- readName}
      }
    #combine two vector containing left flanking and right flanking for single fosmid mapping
    uniq_list_tmp<-c(uniq_list_1,uniq_list_2)
    #make mapped fosmid flankings read name as a unique list, extract read name excluding _L or _R
    uniq_list <-list(unique(substr(uniq_list_tmp,start = 1 , stop = 36)))
    fosmid_counter = fosmid_counter + 1
    names(uniq_list) <- paste0("fos",fosmid_counter)
    fos_list <- c(fos_list,uniq_list)
    #finished make fosmid list in single barcode
  }
  close(readInList)
  fin_fos_list = list()
  R_counter = 1
  #make fosmid list with at least cutoff number of reads
  for (list in fos_list){
    if(length(list) > cutoff){
      shrink_list <- list(list)
      names(shrink_list) <- paste0("R_fos",R_counter)
      fin_fos_list <-c(fin_fos_list , shrink_list)
      R_counter = R_counter + 1
    }else{next}}
  fos_tmp_number = c( paste0("barcode",formatC(j, width = 2, flag = '0')),length(fos_list),length(fin_fos_list))
  fos_sep_df[nrow(fos_sep_df)+1,]= fos_tmp_number
}



#########end of script#################
length(fos_list)
length(fin_fos_list)
View(fos_list)
View(fin_fos_list)
View(fos_sep_df)

#collect information
fos_sep_df = data.frame(barcode = character(), list = integer(), cutoff = integer(),stringsAsFactors=FALSE)
colnames(fos_sep_df) = c("barcode", "list", "list_cutoff")
fos_tmp_number = c( paste0("barcode",formatC(j, width = 2, flag = '0')),length(fos_list),length(fin_fos_list))
fos_sep_df[nrow(fos_sep_df)+1,]= fos_tmp_number

