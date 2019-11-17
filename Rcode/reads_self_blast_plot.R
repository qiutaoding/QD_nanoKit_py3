setwd("F:/nanopore_data/telomere/plot/")
library("seqinr")
library("ggplot2")
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
fileName=paste0("XR/N2L1_1.csv")
df_file=read.table(paste(fileName),header=FALSE)
colnames(df_file)<-blast_title
read1_length=getLength(read.fasta(paste0("XR/N2L1_1.fa")))
read2_length=getLength(read.fasta(paste0("XR/N2L1_1.fa")))
length = 67000
readname = "b0ceb99d-b39b-4c07-8eee-4d31762bbc60"
p<-ggplot()+
  geom_segment(data=df_file,
                mapping=aes(x=qstart,y=sstart,xend=qend,yend=send),
                size=1,
                color="black",linetype=1)+
  xlab(paste(readname))+
  ylab(paste(readname))+
  theme(axis.title=element_text(size=10))+
  scale_x_continuous(breaks=seq(0,length,by=1000),limits=c(0,length),expand=c(0,0),labels=scales::comma)+
  scale_y_continuous(breaks=seq(0,length,by=1000),limits=c(0,length),expand=c(0,0),labels=scales::comma)+
  geom_segment(aes(x=read1_length,y=0,xend=read1_length,yend=read2_length),color="red",linetype="dashed",size=0.5)+
  geom_segment(aes(x=0,y=read2_length,xend=read1_length,yend=read2_length),color="red",linetype="dashed",size=0.5)+
  theme( panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p

pdf(paste0("XR/XR_N2L1_1_plot.pdf"),width=10,height = 10)
print(p)
dev.off()





df_rga =read.table("IIR/rga_N2L1_2.csv",header=FALSE)
colnames(df_rga)<-blast_title
q<- p +
  geom_segment(data=df_rga,
               mapping=aes(x=sstart,y=sstart,xend=send,yend=send),
               size=2,
               color="red",linetype=1)
q  