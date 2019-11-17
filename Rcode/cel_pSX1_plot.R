setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")

df=read.table("psx1/N2EM_blast.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="psx1/N2EM_psx1_read.fa")
readlist <-read.table("psx1/N2EM_psx1_list")
readplot <- function(list_element){
  readName <-as.character(list_element)
  readLength <-length(fa[readName][[1]])
  if (readLength >10000){
    df_tmp <-subset(df,subject == readName)
    df_left <-subset(df_tmp, query =="left" & length >1000)
    df_psx1 <-subset(df_tmp, query =="pSX1" & length >24)
    df_psx1$colour <- sqrt(df_psx1$mismatch+df_psx1$gap+1)
    df_right <- subset(df_tmp, query =="right" &length >100)
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
                        size=1,color=df_psx1$colour,
                        linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))      }
    if(nrow(df_right)>0){
      p<-p+geom_segment(data=df_right,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=1,color="blue", 
                        linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
    if (nrow(df_left)>0 & nrow(df_psx1)> 0 ){
      p<- p+      
        scale_y_continuous(expand = c(0, 0),limits = c(0,2100),breaks=seq(0,2100,by=100))
      pdf(paste0("psx1/plot/N2EM/left/",readName,".pdf"),width=readLength/1000,height = calheight)
      print(p)
      dev.off()  
    }else if(nrow(df_right)>0 & nrow(df_psx1)> 0 ){
      p<- p+      
        scale_y_continuous(expand = c(0, 0),limits = c(0,2100),breaks=seq(0,2100,by=100))
      pdf(paste0("psx1/plot/N2EM/right/",readName,".pdf"),width=readLength/1000,height = calheight)
      print(p)
      dev.off()
    }else if( nrow(df_left)>0 &nrow(df_right)>0  ){
      p<- p+      
        scale_y_continuous(expand = c(0, 0),limits = c(0,2100),breaks=seq(0,2100,by=100))
      pdf(paste0("psx1/plot/N2EM/cross/",readName,".pdf"),width=readLength/1000,height = calheight)
      print(p)
      dev.off()
    }else if( nrow(df_psx1)>20 ){
      p<- p+      
        scale_y_continuous(expand = c(0, 0),limits = c(0,200),breaks=seq(0,200,by=10))
      pdf(paste0("psx1/plot/N2EM/psx1_only/",readName,".pdf"),width=readLength/1000,height = 7)
      print(p)
      dev.off()
    }
  }else{}
}
lapply(readlist$V1, readplot)

