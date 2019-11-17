setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/N2_ChrV_psx1/")
library("ggplot2")
library("seqinr")

df=read.table("N2YA_blast.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="N2YA_V_read.fa")
readlist <-read.table("N2YA_V_list")
readplot <- function(list_element){
  readName <-as.character(list_element)
  readLength <-length(fa[readName][[1]])
  df_tmp <-subset(df,subject == readName)
  df_left <-subset(df_tmp, query =="left" & length >100)
  df_psx1 <-subset(df_tmp, query =="pSX1" & length >24)
  df_right <- subset(df_tmp, query =="right_lgc-26" &length >1000)
  calheight = 15
  p<-ggplot()+
    scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank() )
  if(nrow(df_left)>0){
    p<-p+geom_segment(data=df_left,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="red", 
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_psx1)>0){
    p<-p+geom_segment(data=df_psx1,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="black",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))      }
  if(nrow(df_right)>0){
    p<-p+geom_segment(data=df_right,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="blue", 
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if (nrow(df_left)>0 & nrow(df_psx1)> 0 ){
    p<- p+      
      scale_y_continuous(expand = c(0, 0),limits = c(0,3200),breaks=seq(0,3200,by=100))
    pdf(paste0("plot/N2YA/left/",readName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()  
  }else if(nrow(df_right)>0 & nrow(df_psx1)> 0 ){
    p<- p+      
      scale_y_continuous(expand = c(0, 0),limits = c(0,3200),breaks=seq(0,3200,by=100))
    pdf(paste0("plot/N2YA/right/",readName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()
  }else if( nrow(df_left)>0 &nrow(df_right)>0  ){
    p<- p+      
      scale_y_continuous(expand = c(0, 0),limits = c(0,3200),breaks=seq(0,3200,by=100))
    pdf(paste0("plot/N2YA/cross/",readName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()
  }else if( nrow(df_psx1)> 10 ){
    p<- p+      
      scale_y_continuous(expand = c(0, 0),limits = c(0,3200),breaks=seq(0,3200,by=100))
    pdf(paste0("plot/N2YA/psx1/",readName,".pdf"),width=readLength/1000,height = 7)
    print(p)
    dev.off()
  }
}
lapply(readlist$V1, readplot)

