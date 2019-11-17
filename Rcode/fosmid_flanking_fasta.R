args<-commandArgs(TRUE)
j=as.numeric(args[1])
library(data.table)
fossep_fn = paste0('fos_sep/sample',formatC(j, width = 2, flag = '0'),'.csv')