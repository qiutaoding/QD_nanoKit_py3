setwd("F:/nanopore_data/telomere")
library("ggplot2")
library("seqinr")
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
df=read.table("IR_csv/telo_N2YA.csv",header=FALSE)
colnames(df)<-blast_title
fa<-read.fasta(file="extract/IR/N2YA_IR_read.fa")
readlist <-read.table("extract/IR/N2YA_IR_list")
readplot <-function(list_element){
  readName <-as.character(list_element)
  df_tmp <-subset(df,query == readName)
  df_18s <-subset(df_tmp, subject =="18s" & length >100)
  df_26s <-subset(df_tmp, subject =="26s" &length >100)
  df_ets <- subset(df_tmp , subject == "ets" &length > 100)
  df_left <- subset(df_tmp , subject == "F31C3.6" &length > 100)
  readLength <-length(fa[readName][[1]])
  calheight = 15
  p<-ggplot()+
    scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,5400),breaks=seq(0,5400,by=500))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank())
  if(nrow(df_18s)>0){
    p<-p+geom_segment(data=df_18s,mapping=aes(x=qstart,y=sstart,xend=qend,yend=send),
                      size=1,color="brown", 
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  
  if(nrow(df_26s)>0){
    p<-p+geom_segment(data=df_26s,mapping=aes(x=qstart,y=sstart,xend=qend,yend=send),
                      size=1,color="red",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_ets)>0){
    p<-p+
      geom_segment(data=df_ets,mapping=aes(x=qstart,y=sstart,xend=qend,yend=send),
                   size=4,color="black",
                   linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))+
      geom_segment(data=df_ets,mapping=aes(x=qstart,y=send,xend=qend,yend=sstart),
                   size=4,color="black",linejoin='mitre')
  }
  if(nrow(df_left)>0){
    p<-p+geom_segment(data=df_left,mapping=aes(x=qstart,y=sstart,xend=qend,yend=send),
                      size=2,color="greenyellow",
                      linejoin='mitre')}
  if (nrow(df_left) >0 ){
  } else if (nrow(df_ets)>0){
    if (df_ets[1,]$sstart >df_ets[1,]$send){
      df_check <- df_ets[which.min(df_ets$qstart),]
      if (df_check$sstart >810 & df_check$sstart < 825){
        pdf(paste0("plot/N2YA/telo/R/",readName,".pdf"),width=readLength/1000,height = calheight)
        print(p)
        dev.off()}
    }else{
      df_check <- df_ets[which.max(df_ets$qend),]
      if (df_check$send >810 & df_check$send < 825){
        pdf(paste0("plot/N2YA/telo/F/",readName,".pdf"),width=readLength/1000,height = calheight)
        print(p)
        dev.off()}}}
}
lapply(readlist$V1, readplot)
