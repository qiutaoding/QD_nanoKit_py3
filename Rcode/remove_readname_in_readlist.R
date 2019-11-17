setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/extracted_read")
fh_list=read.table("N2EM_list",header=FALSE,stringsAsFactors = F)
rm_list=read.table("N2EM_5s_only_list",header=FALSE,stringsAsFactors = F)
fh_filter <-fh_list[!(fh_list$V1) %in% rm_list$V1,]

lapply(fh_filter, write, "exclude_5s_only/N2EM_5s_rm_list", append=TRUE, ncolumns=1)
