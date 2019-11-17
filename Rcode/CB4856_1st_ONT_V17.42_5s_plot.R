setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")
df=read.table("plot/CB4856_1st_V17_42.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="extracted_read/CB4856_1st_5s_read.fa")
readlist <-read.table("extracted_read/CB4856_1st_5s_list")
readplot <-function(list_element){
  readName <-as.character(list_element)
  df_tmp <-subset(df,subject == readName)
  df_5s_1 <-subset(df_tmp, query =="cel-5s-unit1" & length >42)
  df_right <- subset(df_tmp , query =="right_boundary" &length > 100)
  df_transition <-subset(df_tmp , query =="transition" &length > 100)
  df_anchor <- subset(df_tmp, query =="T27C5.10" &length > 100)
  df_anchor2 <- subset(df_tmp, query =="srbc-27" &length > 100)
  readLength <-length(fa[readName][[1]])
  calheight = 8
  p<-ggplot()+
    scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,1700),breaks=seq(0,1700,by=100))+
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
  if(nrow(df_transition)>0){
    p<-p+geom_segment(data=df_transition,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=2,color="greenyellow",
                      linejoin='mitre')}
  if(nrow(df_right)>0){
    p<-p+geom_segment(data=df_right,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="blue",
                      linejoin='mitre')}
  if(nrow(df_anchor)>0){
    p<-p+geom_segment(data=df_anchor,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="turquoise",
                      linejoin='mitre')}
  if(nrow(df_anchor2)>0){
    p<-p+geom_segment(data=df_anchor2,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="red",
                      linejoin='mitre')}

  if (nrow(df_anchor) >0 & nrow(df_anchor2)>0 & nrow(df_5s_1)>0){
    pdf(paste0("plot/CB4856_1st/V17_42/",readName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()  
  }else{}
}
lapply(readlist$V1, readplot)

