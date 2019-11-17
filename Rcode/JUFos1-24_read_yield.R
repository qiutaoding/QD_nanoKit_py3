setwd("F:/nanopore_data/JUfosmid")

df=read.table("JUFos1-24.fastq.fai",header=FALSE)
df_short <- df[c("V1","V2")]
library(ggplot2)
bin_set=1000
#make a weighted plot for nanopore sequencing data (data yield at read different length)

p2<-ggplot(df_short,aes(x=V2))+
  geom_histogram(aes(weight = V2),binwidth=1000)+
  scale_x_continuous(labels = scales::comma,expand = c(0, 0),limits = c(-500,90000),breaks=seq(0,90000,10000))+
  scale_y_continuous(labels = scales::comma,expand = c(0, 0),breaks=seq(0,1500000000,250000000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("data")+
  xlab("Read length: binwidth = 1 kb")
p2
pdf("JUFos1-24_data_weight.pdf",width=10,height = 5)
print(p2)
dev.off()
