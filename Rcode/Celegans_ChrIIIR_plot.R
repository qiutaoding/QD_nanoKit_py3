setwd("F:/nanopore_data/telomere/")
library("ggplot2")
library("seqinr")
df=read.table("blastn_check/N2YA_IIIR.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="extract/IIIR/history/N2YA_IIIR_read.fa")
readlist <-read.table("extract/IIIR/history/N2YA_IIIR_list")
readplot <-function(list_element){
  readName <-as.character(list_element)
  df_tmp <-subset(df,subject == readName)
  df_left <-subset(df_tmp, query =="pdr1" & length >100)
  df_mdt <-subset(df_tmp, query =="mdt-29" & length >100)
  df_pot3 <-subset(df_tmp, query =="pot3" & length >100)
  readLength <-length(fa[readName][[1]])
  calheight = 10
  fileName <-substring(readName,19)
  p<-ggplot()+
    scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
    scale_y_continuous(expand = c(0, 0),limits = c(0,5200),breaks=seq(0,5200,by=200))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  if(nrow(df_left)>0){
    p<-p+geom_segment(data=df_left,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="black", 
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}

  if(nrow(df_mdt)>0){
    p<-p+geom_segment(data=df_mdt,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="blue",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_pot3)>0){
    p<-p+geom_segment(data=df_pot3,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=2,color="red",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}


  if (nrow(df_left) >0 & (nrow(df_mdt)==0 & nrow(df_pot3)==0)){
    pdf(paste0("plot/IIIR/left_only/",fileName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()  
  }else if (nrow(df_left) >0){
    pdf(paste0("plot/IIIR/left/",fileName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()
  }else {
    pdf(paste0("plot/IIIR/un/",fileName,".pdf"),width=readLength/1000,height = calheight)
    print(p)
    dev.off()
  } 
}
lapply(readlist$V1, readplot)

