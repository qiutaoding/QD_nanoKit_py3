setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/N2_18S/blast/")
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")

gC_EM <- 100.272
gC_L1 <-35.6229
gC_YA <-48.2667

df_EM=read.table("N2EM_5_8s.csv",header=FALSE)
df_L1=read.table("N2L1_5_8s.csv",header=FALSE)
df_YA=read.table("N2YA_5_8s.csv",header=FALSE)
colnames(df_EM)<-blast_title
colnames(df_L1)<-blast_title
colnames(df_YA)<-blast_title
cutoff <- 22
gene_length <-153
filtter_EM <- subset(df_EM,length>cutoff)
filtter_L1 <- subset(df_L1,length>cutoff)
filtter_YA <- subset(df_YA,length>cutoff)

sum(abs(filtter_EM$send-filtter_EM$sstart)+1)/gene_length/gC_EM
sum(abs(filtter_L1$send-filtter_L1$sstart)+1)/gene_length/gC_L1
sum(abs(filtter_YA$send-filtter_YA$sstart)+1)/gene_length/gC_YA


