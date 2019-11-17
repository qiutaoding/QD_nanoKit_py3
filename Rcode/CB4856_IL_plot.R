setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")
df=read.table("CB4856_extract/CB4856_IL_psx1.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="CB4856_extract/CB4856_IL_read.fa")
readlist <-read.table("CB4856_extract/CB4856_IL_list")
readplot <- function(list_element){
  readName <-as.character(list_element)
  readLength <-length(fa[readName][[1]])
  if (readLength >10000){
    df_tmp <-subset(df,subject == readName)
    df_unit <-subset(df_tmp, query =="n2_rdna_unit" & length >100)
    df_26s <-subset(df_tmp, query =="26s" & length >100)
    df_psx1 <-subset(df_tmp, query =="pSX1" & length >24)
    df_fragment <- subset(df_tmp, query =="f6c1963e-b5ef-4d14-8b45-2d2aa61d6ea3:48076-59645" &length >100)
    calheight = 15
    p<-ggplot()+
      scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
      scale_y_continuous(expand = c(0, 0),limits = c(0,12000),breaks=seq(0,12000,by=1000))+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank() )
    if(nrow(df_unit)>0){
      p<-p+geom_segment(data=df_unit,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=1,color="brown", 
                        linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
    if(nrow(df_psx1)>0){
      p<-p+geom_segment(data=df_psx1,mapping=aes(x=sstart,y=qstart*10,xend=send,yend=qend*10),
                        size=1,color="black",
                        linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
    if(nrow(df_fragment)>0){
      p<-p+geom_segment(data=df_fragment,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=1,color="red",
                        linejoin='mitre')}
    if(nrow(df_26s)>0){
      p<-p+geom_segment(data=df_26s,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                        size=2,color="blue", 
                        linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}

    
    if (nrow(df_unit)>2 & nrow(df_fragment) >0 & nrow(df_psx1) >10){
      pdf(paste0("CB4856_extract/psx1/",readName,".pdf"),width=readLength/1000,height = calheight)
      print(p)
      dev.off()  
    }else{}
  }else{}
}
lapply(readlist$V1, readplot)

