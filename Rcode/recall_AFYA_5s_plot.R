setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")
df=read.table("plot/AFYA_5s.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title
fa<-read.fasta(file="extracted_read/AFYA_read.fa")
readlist <-read.table("extracted_read/AFYA_list")
readplot <-function(list_element){
  readName <-as.character(list_element)
  df_tmp <-subset(df,subject == readName)
  df_5s_1 <-subset(df_tmp, query =="cbr-unit-1" & length >200)
  df_5s_2 <-subset(df_tmp, query =="cbr-unit-2" & length >200)
  df_5s_3 <-subset(df_tmp, query =="cbr-unit-2-short" & length >500)
  df_CBG06809 <- subset(df_tmp, query =="CBG06809" &length > 100)
  df_CBG06809_gap <- subset(df_tmp , query == "CBG06809_gap" &length > 100)
  df_transition <-subset(df_tmp , query =="transition_gap" &length > 100)
  df_right <- subset(df_tmp , query =="CBG10685" &length > 100)
  readLength <-length(fa[readName][[1]])
  calheight = 8
  p<-ggplot()+
    xlab(paste(readName,'length:',readLength,'bp'))+
    scale_x_continuous(expand = c(0, 0),breaks=seq(0,readLength,by=1000),labels=scales::comma,limits = c(0,readLength))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  if(nrow(df_5s_1)>0){
    p<-p+geom_segment(data=df_5s_1,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=2,color="black", 
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_5s_2)>0){
    p<-p+geom_segment(data=df_5s_2,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="red",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_5s_3)>0){
    p<-p+geom_segment(data=df_5s_3,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="navy",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_CBG06809)>0){
    p<-p+geom_segment(data=df_CBG06809,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=2,color="red",
                      linejoin='mitre')}
  if(nrow(df_CBG06809_gap)>0){
    p<-p+geom_segment(data=df_CBG06809_gap,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=2,color="greenyellow",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_right)>0){
    p<-p+geom_segment(data=df_right,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="blue",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if(nrow(df_transition)>0){
    p<-p+geom_segment(data=df_transition,mapping=aes(x=sstart,y=qstart,xend=send,yend=qend),
                      size=1,color="turquoise",
                      linejoin='mitre',arrow = arrow(length = unit(0.5, "cm")))}
  if ((nrow(df_CBG06809) >0 |nrow(df_CBG06809_gap) >0 )& nrow(df_5s_1) ==0){
    p<- p+scale_y_continuous(expand = c(0, 0),limits = c(0,2600),breaks=seq(0,2600,by=100))
    pdf(paste0("plot/AFYA/left_only/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off()  
  }else if ((nrow(df_CBG06809) >0 |nrow(df_CBG06809_gap) > 0 ) & (nrow(df_5s_1) > 0|nrow(df_5s_2) > 0)){
    p<-p+scale_y_continuous(expand = c(0, 0),limits = c(0,2600),breaks=seq(0,2600,by=100))
    pdf(paste0("plot/AFYA/left/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off()  
  }else if (nrow(df_right)>0 &(nrow(df_5s_1)>0|nrow(df_5s_2) > 0 |nrow(df_5s_3) > 0)){
    p<- p+scale_y_continuous(expand = c(0, 0),limits = c(0,1400),breaks=seq(0,1400,by=100))
    pdf(paste0("plot/AFYA/right/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off() 
  }else if (nrow(df_right)>0 &(nrow(df_5s_1)==0|nrow(df_5s_2) == 0|nrow(df_5s_3) == 0)&nrow(df_transition)==0){
    p<- p+scale_y_continuous(expand = c(0, 0),limits = c(0,1400),breaks=seq(0,1400,by=100))
    pdf(paste0("plot/AFYA/right_only/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off() 
  }else if (nrow(df_transition)>0 &(nrow(df_5s_1)>0|nrow(df_5s_2) > 0|nrow(df_5s_3) > 0)){
    p<- p+scale_y_continuous(expand = c(0, 0),limits = c(0,1400),breaks=seq(0,1400,by=100))
    pdf(paste0("plot/AFYA/transition/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off() 
  }else if (nrow(df_5s_1)>0 & nrow(df_5s_2) > 0){
    p<-p+scale_y_continuous(expand = c(0, 0),limits = c(0,1000),breaks=seq(0,1000,by=100))
    pdf(paste0("plot/AFYA/combine_5s_1_2/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off() 
  }else if (nrow(df_5s_1)>0 & (nrow(df_5s_2) == 0|nrow(df_5s_3) == 0)){
    p<- p+scale_y_continuous(expand = c(0, 0),limits = c(0,1000),breaks=seq(0,1000,by=100))
    pdf(paste0("plot/AFYA/5s_1_only/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off() 
  }else if (nrow(df_5s_1)==0 & (nrow(df_5s_2) > 0|nrow(df_5s_3) > 0)){
    p<- p+scale_y_continuous(expand = c(0, 0),limits = c(0,1000),breaks=seq(0,1000,by=100))
    pdf(paste0("plot/AFYA/5s_2/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off() 
  }else{
    pdf(paste0("plot/AFYA/un/",readName,".pdf"),width=5+readLength/1000,height = calheight)
    print(p)
    dev.off() 
  }
}
lapply(readlist$V1, readplot)
