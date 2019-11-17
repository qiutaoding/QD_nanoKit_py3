setwd("F:/nanopore_data/telomere/extract/CB4856/check_dup/")
library("seqinr")
library("ggplot2")
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
location = "IIL_3"
df_file=read.table(paste0(location,".csv"),header=FALSE)
colnames(df_file)<-blast_title
read1_length=getLength(read.fasta(paste0(location,".fa")))
read2_length=getLength(read.fasta(paste0(location,".fa")))
length = (as.integer(read1_length/1000)+1)*1000
readname = "7b4d48cf-883e-4c84-96eb-60aa5963a892"
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
  #rDNA region
  #geom_segment(aes(x=0,y=17154,xend=17154,yend=17154),color="blue",size=0.5)+
  #geom_segment(aes(x=17154,y=0,xend=17154,yend=17154),color="blue",size=0.5)+
  #Telomere region
  #geom_segment(aes(x=96571,y=96571,xend=96571,yend=read1_length),color="blue",size=0.5)+
  #geom_segment(aes(x=96571,y=96571,xend=read1_length,yend=96571),color="blue",size=0.5)+
  #annotation
  #annotate("text", x= 7000, y= 7000, label= "45S rDNA",color = "red", size = 8)+
  #annotate("text", x= 100000, y= 100000, label= "Telomere",color = "white", size = 4)+
  theme( panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p

pdf(paste0("IIL_3_plot.pdf"),width=10,height = 10)
print(p)
dev.off()
