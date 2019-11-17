setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")
df=read.table("plot/AFYA_5s_20190811.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
#read aligned length
df$qlength <-df$qend-df$qstart+1
#extract unique aligned
df_1 <-subset(df,subject =="cbr-unit-1" & qlength>170)
df_2 <-subset(df,subject =="cbr-unit-2" &qlength>170)
gcov=90.8984
unit1_CN <- sum(df_1$qlength)/938/gcov
unit2_CN <- sum(df_2$qlength)/693/gcov
unit1_CN
unit2_CN

al_1_CN = sum(df_1$length)/938/gcov
al_2_CN = sum(df_2$length)/693/gcov
al_1_CN
al_2_CN


r_1_CN = sum(abs(df_1$sstart-df_1$send))/938/gcov
r_2_CN = sum(abs(df_2$sstart-df_2$send))/693/gcov
r_1_CN
r_2_CN
