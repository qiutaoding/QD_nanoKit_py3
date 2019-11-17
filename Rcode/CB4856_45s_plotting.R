setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")
df=read.table("plot/CB4856_45s.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="extracted_read/CB4856_45s_read.fa")
readlist <-read.table("extracted_read/CB4856_45s_list")
readplot <- function(list_element){
  readName <-as.character(list_element)
  readLength <-length(fa[readName][[1]])
  if (readLength >10000){
    df_tmp <-subset(df,subject == readName)
    df_18s <-subset(df_tmp, query =="18s" & length >100)
    df_5_8s <- subset(df_tmp, query =="5_8s" & length >40)
    df_26s <-subset(df_tmp, query =="26s" &length >100)
    df_ets <- subset(df_tmp , query == "ets" &length > 100)
    df_left <- subset(df_tmp , query == "F31C3.6" &length > 500)
    df_telo <- subset(df_tmp , query == "F31C3.6" & qstart >3762 & qend <3887 )
    df_unit <- subset(df_tmp , query == "rnda_unit" &length > 100)
    calheight = 15
    p<-ggplot()+
      scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
      scale_y_continuous(expand = c(0, 0),limits = c(0,7200),breaks=seq(0,7200,by=100))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank() )
    if(nrow(df_18s)>0){
      p<-p+geom_segment(data=df_18s,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=1,color="black", 
                        linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
    if(nrow(df_26s)>0){
      p<-p+geom_segment(data=df_26s,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=1,color="red",
                        linejoin='mitre')}
    if(nrow(df_ets)>0){
      p<-p+geom_segment(data=df_ets,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=1,color="green",linejoin='mitre')}
    if(nrow(df_left)>0){
      p<-p+geom_segment(data=df_left,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=2,color="blue",linejoin='mitre')}
    if(nrow(df_telo)>0){
      p<-p+geom_segment(data=df_telo,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=2,color="blue",linejoin='mitre')}
    if(nrow(df_unit)>0){
      p<-p+geom_segment(data=df_unit,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=4,color="black",linejoin='mitre')}
    
    if (nrow(df_left)>0 & nrow(df_26s) >0 ){
      pdf(paste0("plot/CB4856_ONT/45S/left/",readName,".pdf"),width=readLength/1000,height = calheight)
      print(p)
      dev.off()  
    }else if (nrow(df_telo)>10 & nrow(df_18s)>2){
      pdf(paste0("plot/CB4856_ONT/45S/telo/",readName,".pdf"),width=readLength/1000,height = calheight)
      print(p)
      dev.off()  
    }else if (nrow(df_18s)>3){
      pdf(paste0("plot/CB4856_ONT/45S/long_45s/",readName,".pdf"),width=readLength/1000,height = calheight)
      print(p)
      dev.off()  
    }else{}
  }else{}
}
lapply(readlist$V1, readplot)

