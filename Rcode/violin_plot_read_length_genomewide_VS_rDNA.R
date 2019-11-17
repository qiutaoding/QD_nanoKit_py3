setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/genomewide_readlength_distribution/")
library("ggplot2")
#####by average length in 999 positions crossing genome####
df=read.table("position_average/cel_genomewide_distribution.csv",header=TRUE,sep=",",stringsAsFactors = F)
p <-ggplot(df,aes(x=stage,y=avg_len,fill = stage))+
  geom_violin()+
  geom_boxplot(width=0.2, notch=TRUE, outlier.size=1)+
  scale_y_continuous(expand = c(0, 0),breaks=seq(0,45000,by=5000),labels= seq(0,45,by=5),limits = c(0,45000))+
  #amanually add 5S rDNA read length average
  geom_point(aes(x="N2EM", y = 14872.5),color = "red", shape= 15, size=3)+
  geom_point(aes(x="N2L1", y = 17290.6),color = "red", shape= 15, size=3)+
  geom_point(aes(x="N2YA", y = 14260.7),color = "red", shape= 15, size=3)+
  #amanually add 45S rDNA read length average
  geom_point(aes(x="N2EM", y = 13975),color = "blue", shape= 17,size=3)+
  geom_point(aes(x="N2L1", y = 17818.3),color = "blue", shape= 17,size=3)+
  geom_point(aes(x="N2YA", y = 12813.9),color = "blue", shape= 17,size=3)+
  scale_x_discrete(labels=c("EMB","L1","YA"))+
  labs(x="Stage", y= "Average read length  (kb)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(size=12, colour = "black"),
        axis.title = element_text(size=14, colour = "black"))
p

#export from rihgt panel is better
pdf("position_average/rDNA_VS_genomewide_read_length_distribution_20191101.pdf",width=4,height = 4)
print(p)
dev.off()

####by individual read length####
library("ggplot2")
df_EM=read.table("N2EM_dis.csv",header=F,sep="\t",stringsAsFactors = F)
colnames(df_EM)<-c("readname", "length","offset","linebase","linewidth")
df_EM$stage <- "N2EM"
df_L1=read.table("N2L1_dis.csv",header=F,sep="\t",stringsAsFactors = F)
colnames(df_L1)<-c("readname", "length","offset","linebase","linewidth")
df_L1$stage <- "N2L1"
df_YA=read.table("N2YA_dis.csv",header=F,sep="\t",stringsAsFactors = F)
colnames(df_YA)<-c("readname", "length","offset","linebase","linewidth")
df_YA$stage <- "N2YA"

df <-rbind(df_EM,df_L1,df_YA)


p <-ggplot(df,aes(x=stage,y=length,fill = stage))+
  geom_violin()+
  geom_boxplot(width=0.1)

p
