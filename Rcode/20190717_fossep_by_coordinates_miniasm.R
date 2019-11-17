args<-commandArgs(TRUE)
j=as.numeric(args[1])
library(dplyr)
library(data.table)
header = c('qname', 'qlen','qstart','qend','strand','rname','rlen','rstart','rend','match','alignlength','mapq',"primary")
trim_info_fn = paste0('align/miniasm_align/trimmed_barcode',formatC(j, width = 2, flag = '0'),'.paf')
trim_info <- read.csv(trim_info_fn,sep = "\t", header = F)
#cutoff session
cutoff = 500 #flanking length cutoff
read_count_cutoff = 1 #minimal reads required in single fosmid location
boundary_cutoff = 300
df <-subset(trim_info, trim_info$V3<boundary_cutoff & (trim_info$V2-trim_info$V4)<boundary_cutoff &trim_info$V12 >= 10 &trim_info$V13 =="tp:A:P",select = c(seq(1,13)) )
colnames(df)<-header
df_sorted <-df[order(df$rname,df$rend),]
df_L <-df_sorted[df_sorted$qname %like% "_L",]
df_R <-df_sorted[df_sorted$qname %like% "_R",]
rm(df)
rm(df_sorted)
fosmid_coordinate = data.frame(fosmid = character(),
                               L_contig = character(),
                               L_start = integer(),
                               L_strand =character(),
                               L_confidence = integer(),
                               R_contig = character(),
                               R_end = integer(),
                               R_strand= character(), 
                               R_confidence= integer(),
                               stringsAsFactors=FALSE)
counter = 0
left_flanking_working_list <- list()
right_flanking_remove_list <- list()
for (read_list in df_L$qname){
  if (read_list %in% left_flanking_working_list){next}else{
    readname = substr(read_list,start =1 ,stop =36)
    tmp_L_df <-df_L[df_L$qname==read_list,]
    tmp_L_df<-tmp_L_df[which.max(tmp_L_df$mapq),]
    tmp_L_df<-tmp_L_df[which.max(tmp_L_df$match),]
    find_similar_R_tmp <-df_R[FALSE,]
    if (tmp_L_df$strand == "+"){
      tmp_L_start = tmp_L_df$rend
      L_contig_info = as.character(tmp_L_df$rname)
      find_similar <- subset(df_L,(strand=="+" &
                                     rname == L_contig_info & 
                                     rend <tmp_L_start +cutoff & 
                                     rend > tmp_L_start-cutoff))
      find_similar <-find_similar[!(find_similar$qname %in% left_flanking_working_list),]
      L_confidence = nrow(find_similar)
      if (nrow(find_similar) < read_count_cutoff){
        next
      }else{
        counter = counter +1
        L_strand_info <- "+"
        L_start_info <-median(find_similar$rend)
        similar_list <-unlist(substr(find_similar$qname,start = 1, stop = 36))
        left_flanking_working_list <-c(left_flanking_working_list,as.character(find_similar$qname))
        for (listname in similar_list){
          if (nrow(df_R[df_R$qname %like% listname, ])!=0){
            df_R_adding <-df_R[df_R$qname %like% listname, ]
            find_similar_R_tmp[nrow(find_similar_R_tmp)+1,] <- df_R_adding[which.max(df_R_adding$mapq),]
          }else{next}}
        if (nrow(find_similar_R_tmp) != 0){        
          R_contig_df <-data.frame(table(find_similar_R_tmp$rname))
          R_strand_df <-data.frame(table(find_similar_R_tmp$strand))
          R_contig_info <-as.character(R_contig_df[which.max(R_contig_df$Freq),]$Var1)
          R_strand_info <-as.character(R_strand_df[which.max(R_strand_df$Freq),]$Var1)
          if (R_strand_info=="+"){
            tmp_R_stop = median(find_similar_R_tmp$rstart)
            find_similar_R <- subset(df_R, rname==R_contig_info & rstart <tmp_R_stop+cutoff & rstart >tmp_R_stop-cutoff)
            right_flanking_remove_list <- c(right_flanking_remove_list, as.character(find_similar_R$qname))
            R_stop_info <-median(find_similar_R$rstart)
            R_confidence =as.integer(R_strand_df[which.max(R_strand_df$Freq),]$Freq)
          }else{
            tmp_R_stop = median(find_similar_R_tmp$rend)
            find_similar_R <- subset(df_R, rname==R_contig_info & rend <tmp_R_stop+cutoff & rend >tmp_R_stop-cutoff)
            right_flanking_remove_list <- c(right_flanking_remove_list, as.character(find_similar_R$qname))
            R_stop_info <-median(find_similar_R$rend)
            R_confidence =as.integer(R_strand_df[which.max(R_strand_df$Freq),]$Freq)}
        }else{
          R_contig_info = "un"
          R_stop_info =0
          R_strand_info = "un"
          R_confidence = 0}
        fosmid_coordinate[nrow(fosmid_coordinate)+1,]<-c(counter,
                                                         L_contig_info,
                                                         as.integer(L_start_info),
                                                         L_strand_info,
                                                         as.integer(L_confidence),
                                                         R_contig_info,
                                                         as.integer(R_stop_info),
                                                         R_strand_info,
                                                         as.integer(R_confidence))
      }
    }else if (tmp_L_df$strand == "-"){
      # if left flanking aligned to negative strand
      tmp_L_start = tmp_L_df$rstart
      L_contig_info = as.character(tmp_L_df$rname)
      find_similar <- subset(df_L,(strand=="-" & 
                                     rname == L_contig_info &
                                     rstart < tmp_L_start +cutoff &
                                     rstart > tmp_L_start -cutoff ))
      find_similar <-find_similar[!(find_similar$qname %in% left_flanking_working_list),]
      L_confidence = nrow(find_similar)
      if (nrow(find_similar) < read_count_cutoff){
        next
      }else{
        counter = counter +1
        L_strand_info <- "-"
        L_start_info <-median(find_similar$rstart)
        similar_list <-unlist(substr(find_similar$qname,start = 1, stop = 36))
        left_flanking_working_list <-c(left_flanking_working_list,as.character(find_similar$qname))
        for (listname in similar_list){
          if (nrow(df_R[df_R$qname %like% listname, ])!=0){
            df_R_adding <-df_R[df_R$qname %like% listname, ]
            find_similar_R_tmp[nrow(find_similar_R_tmp)+1,] <- df_R_adding[which.max(df_R_adding$mapq),]
          }else{next}}
        if (nrow(find_similar_R_tmp) != 0){
          R_contig_df <-data.frame(table(find_similar_R_tmp$rname))
          R_strand_df <-data.frame(table(find_similar_R_tmp$strand))
          R_contig_info <-as.character(R_contig_df[which.max(R_contig_df$Freq),]$Var1)
          R_strand_info <-as.character(R_strand_df[which.max(R_strand_df$Freq),]$Var1)
          if (R_strand_info=="+"){
            tmp_R_stop = median(find_similar_R_tmp$rstart)
            find_similar_R <- subset(df_R, rname==R_contig_info & rstart <tmp_R_stop+cutoff & rstart >tmp_R_stop-cutoff)
            right_flanking_remove_list <- c(right_flanking_remove_list, as.character(find_similar_R$qname))
            R_stop_info <-median(find_similar_R$rstart)
            R_confidence =as.integer(R_strand_df[which.max(R_strand_df$Freq),]$Freq)
          }else{
            tmp_R_stop = median(find_similar_R_tmp$rend)
            find_similar_R <- subset(df_R, rname==R_contig_info & rend <tmp_R_stop+cutoff & rend >tmp_R_stop-cutoff)
            right_flanking_remove_list <- c(right_flanking_remove_list, as.character(find_similar_R$qname))
            R_stop_info <-median(find_similar_R$rend)
            R_confidence =as.integer(R_strand_df[which.max(R_strand_df$Freq),]$Freq)}
        }else{
          R_contig_info = "un"
          R_stop_info =0
          R_strand_info = "un"
          R_confidence = 0}
        fosmid_coordinate[nrow(fosmid_coordinate)+1,]<-c(counter,
                                                         L_contig_info,
                                                         as.integer(L_start_info),
                                                         L_strand_info,
                                                         as.integer(L_confidence),
                                                         R_contig_info,
                                                         as.integer(R_stop_info),
                                                         R_strand_info,
                                                         as.integer(R_confidence))
      }
    }
  }
}  
df_R_remain <- df_R[!(df_R$qname %in% right_flanking_remove_list),]
right_flanking_working_list <-list()
for (read_list in df_R_remain$qname){
  tmp_R_df <-df_R_remain[df_R_remain$qname==read_list,]
  tmp_R_df<- tmp_R_df[which.max(tmp_R_df$mapq),]
  tmp_R_df<-tmp_R_df[which.max(tmp_R_df$match),]
  if (tmp_R_df$strand == "+"){
    tmp_R_start = tmp_R_df$rstart
    R_contig_info = as.character(tmp_R_df$rname)
    find_similar <- subset(df_R_remain,(strand=="+" & rname == R_contig_info & rstart <tmp_R_start+cutoff & rstart > tmp_R_start-cutoff))
    find_similar <-find_similar[!(find_similar$qname %in% right_flanking_working_list),]
    R_confidence <- nrow(find_similar)
    if (nrow(find_similar) < read_count_cutoff){
      next
    }else{
      counter = counter +1
      R_strand_info <- "+"
      R_stop_info <-median(find_similar$rstart)
      similar_list <-unlist(substr(find_similar$qname,start = 1, stop = 36))
      right_flanking_working_list <-c(right_flanking_working_list,as.character(find_similar$qname))
      L_contig_info = "un"
      L_start_info =0
      L_strand_info = "un"
      L_confidence =0
      fosmid_coordinate[nrow(fosmid_coordinate)+1,]<-c(counter,
                                                       L_contig_info,
                                                       as.integer(L_start_info),
                                                       L_strand_info,
                                                       as.integer(L_confidence),
                                                       R_contig_info,
                                                       as.integer(R_stop_info),
                                                       R_strand_info,
                                                       as.integer(R_confidence))
    }
  }else if (tmp_R_df$strand == "-"){
    tmp_R_start = tmp_R_df$rend
    R_contig_info = as.character(tmp_R_df$rname)
    find_similar <- subset(tmp_R_df,(strand=="-" & rname == R_contig_info & rend <tmp_R_start+cutoff & rend > tmp_R_start-cutoff))
    find_similar <-find_similar[!(find_similar$qname %in% right_flanking_working_list),]
    R_confidence <- nrow(find_similar)
    if (nrow(find_similar) < read_count_cutoff){
      next
    }else{
      counter = counter +1
      R_strand_info <- "-"
      R_stop_info <-median(find_similar$rend)
      similar_list <-unlist(substr(find_similar$qname,start = 1, stop = 36))
      left_flanking_working_list <-c(right_flanking_working_list,as.character(find_similar$qname))
      L_contig_info = "un"
      L_start_info =0
      L_strand_info = "un"
      L_confidence =0
      fosmid_coordinate[nrow(fosmid_coordinate)+1,]<-c(counter,
                                                       L_contig_info,
                                                       as.integer(L_start_info),
                                                       L_strand_info,
                                                       as.integer(L_confidence),
                                                       R_contig_info,
                                                       as.integer(R_stop_info),
                                                       R_strand_info,
                                                       as.integer(R_confidence))
    }
  }
} 

write.table(data.frame(fosmid_coordinate),
            paste0('fos_by_coordinate/miniasm/barcode',formatC(j,width=2,flag='0'),'.csv'),
            col.names=T,row.names = F,sep='\t')
cat(paste0("Complete barcode",formatC(j,width=2,flag='0')," and write ",counter, " fosmid into file.\n"))
