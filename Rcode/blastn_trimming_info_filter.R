args<-commandArgs(TRUE)
j=as.numeric(args[1])
header = c('qname', 'rname','percent','length','mis','gap','qstart','qend','rstart','rend','evalue','bit')
trim_info_fn = paste0('blastn/barcode',formatC(j, width = 2, flag = '0'),'.csv')
trim_info <- read.csv(trim_info_fn,sep = "\t", header = F)
colnames(trim_info)<-header
read_list_fn = paste0('blastn/barcode',formatC(j, width = 2, flag = '0'),'.list')
read_list <- read.csv(read_list_fn,sep = "\t", header = F)
read_trim_df = data.frame(readname = character(), start=integer(),end =integer(),strand=character(), stringsAsFactors=FALSE)
saving = paste0('blastn_backbone_info/barcode',formatC(j, width = 2, flag = '0'),'.csv')
merge_blast <-function(read_name){
  readname=as.character(read_name)
  tmp<-subset(trim_info, qname==readname)
  start = min(tmp$qstart)
  end = max(tmp$qend)
  strand_info = tmp[which.max(tmp$length),]
  if (strand_info$rend-strand_info$rstart >0){
    strand = "+"
  }else{strand ="-"}
  data.frame(readname,start,end,strand)
}
read_trim_df <- lapply(read_list$V1,merge_blast)
lapply(1:length(read_trim_df), function(x) write.table(read_trim_df[[x]],
                                                       file=saving,
                                                       sep='\t',
                                                       row.names=FALSE,
                                                       col.names = FALSE, 
                                                       append=T))

print(paste0('==========Finished barcode',formatC(j, width = 2, flag = '0'),'=========='))

