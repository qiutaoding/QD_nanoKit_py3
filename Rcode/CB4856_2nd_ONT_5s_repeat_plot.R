setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")
df=read.table("CB4856_extract/repeat_blastn.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="extracted_read/CB4856_2nd_5s_read.fa")
readlist <-read.table("extracted_read/CB4856_2nd_5s_list")
readplot <- function(list_element){
  readName <-as.character(list_element)
  df_tmp <-subset(df,subject == readName)
  df_5s_1 <-subset(df_tmp, query =="cel-5s-unit1" & length >20)
  df_repeat <-subset(df_tmp, query =="repeat_from_N2_genome" & length >100)
  df_T27C5_10 <- subset(df_tmp, query =="T27C5.10" &length > 100)
  readLength <-length(fa[readName][[1]])
  calheight = 15
  p<-ggplot()+
    scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,7000),breaks=seq(0,7000,by=100))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank() )
  if(nrow(df_5s_1)>0){
    p<-p+geom_segment(data=df_5s_1,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="black", 
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_repeat)>0){
    p<-p+geom_segment(data=df_repeat,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=2,color="red",
                      linejoin='mitre')}

  if(nrow(df_T27C5_10)>0){
    p<-p+geom_segment(data=df_T27C5_10,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="violetred1",linejoin='mitre')}
  
  if (nrow(df_repeat)>1 & nrow(df_5s_1) >4 ){
    pdf(paste0("plot/CB4856_repeat/left/",readName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()  
  }else if (nrow(df_T27C5_10) >0 & nrow(df_5s_1)>0 &nrow(df_repeat)>0){
    pdf(paste0("plot/CB4856_repeat/V17_42/",readName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()  
  }else if (nrow(df_repeat)>0 & nrow(df_5s_1) >0){
    pdf(paste0("plot/CB4856_repeat/right/",readName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()
  }else{}
}
lapply(readlist$V1, readplot)

