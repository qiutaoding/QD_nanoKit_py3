setwd("F:/nanopore_data/rDNA/")
library("ggplot2")
df=read.table("N2_YA_gDNA.fastq.fai",header=FALSE)
colnames(df)<-c("readname", "length","offset","linebase","linewidth","qual")
#statistics
summary(df$length)

#read length distribution
p <-ggplot(df,aes(x=length))+ 
  geom_histogram(binwidth=1000)+
  scale_x_continuous(expand = c(0, 0),breaks=seq(0,70000,by=10000),limits = c(0,70000),labels=scales::comma)+
  #scale_y_continuous(trans='log2',expand = c(0, 0),breaks=c(2**(1:14)))+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
p
pdf("recall_rdna_analysis/5S_flanking_depth/N2YA_reads_distribution.pdf",width=10,height = 4)
print(p)
dev.off()




df=read.table("zzy0603_gDNA_20190914.fastq.fai",header=FALSE)
colnames(df)<-c("readname", "length","offset","linebase","linewidth","qual")
#read length distribution
p <-ggplot(df,aes(x=length))+ 
  geom_histogram(binwidth=1000)+
  scale_x_continuous(expand = c(0, 0),breaks=seq(0,100000,by=10000),limits = c(0,100000),labels=scales::comma)+
  #scale_y_continuous(trans='log2',expand = c(0, 0),breaks=c(2**(1:14)))+
  scale_y_continuous(expand = c(0, 0),labels=scales::comma)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
p
pdf("recall_rdna_analysis/5S_flanking_depth/zzy0603_reads_distribution.pdf",width=10,height = 4)
print(p)
dev.off()

