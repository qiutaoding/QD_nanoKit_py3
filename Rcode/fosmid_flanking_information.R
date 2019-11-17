args<-commandArgs(TRUE)
j=as.numeric(args[1])
library(data.table)
fossep_fn = paste0('fos_sep/sample',formatC(j, width = 2, flag = '0'),'.csv')
fos_sep <- read.csv(fossep_fn,sep = "\t", header = T)
fasta_fh <-read.csv(paste0("trimmed/trimmed_barcode",formatC(j, width = 2, flag = '0'),".fa.fai"),sep = "\t",header =F)
fasta_index <- subset(fasta_fh, select = c(1,2))
colnames(fasta_index) <-c("readname","length")
fos_sep$L_length <-NA
fos_sep$R_length <-NA
rm(fasta_fh)
fos_flanking_info_df = data.frame(fos_number = integer(),
                                  left_flanking=character(),
                                  right_flanking=character(),  
                                  stringsAsFactors=FALSE)

for (i in seq(1,max(fos_sep$fos_number))){
  tmp_df <-subset(fos_sep, fos_number == i)
  for (read_name in tmp_df$read){
    if (nrow(fasta_index[fasta_index$readname %like% paste0(read_name,"_L"),]) == 1){
      tmp_df[tmp_df$read == read_name,]$L_length <-fasta_index[fasta_index$readname %like% paste0(read_name,"_L"),]$length
    }
    if (nrow(fasta_index[fasta_index$readname %like% paste0(read_name,"_R"),]) == 1){
      tmp_df[tmp_df$read == read_name,]$R_length <-fasta_index[fasta_index$readname %like% paste0(read_name,"_R"),]$length
    }
  }
  left_read = as.character(tmp_df[which.max(tmp_df$L_length),]$read)
  right_read = as.character(tmp_df[which.max(tmp_df$R_length),]$read)
  if (length(left_read) == 0 ){left_read = "empty"}
  if (length(right_read) == 0 ){right_read = "empty"}
  fos_flanking_info_df[nrow(fos_flanking_info_df)+1,]<-c(i,left_read,right_read )
}
write.table(data.frame(fos_flanking_info_df),
            paste0('fosmid_flanking_infocsv/sample',formatC(j,width=2,flag='0'),'.csv'),
            col.names=T,row.names = F,sep='\t')


