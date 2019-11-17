setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/")
library("ggplot2")
library("seqinr")
df=read.table("plot/N2YA_5s.csv",header=FALSE)
blast_title<-c("query","subject","ldentity","length","mismatch","gap","qstart","qend","sstart","send","evalue","bitscore")
colnames(df)<-blast_title

##################N2YA_left######################
left_list <- read.table("plot/N2YA/left_read.list")
left_df <- data.frame(readID = character(), 
                      remain= integer(),
                      strand = character(), 
                      stringsAsFactors = FALSE)
for (left_read in left_list$V1){
  readName <-as.character(left_read)
  df_tmp <-subset(df,df$subject == readName)
  df_5s_1 <-subset(df_tmp, query =="cel-5s-unit1" & length >42)
  df_1st_gap <- subset(df_tmp , query == "First_gap" &length > 100)
  df_left_8 <- subset(df_tmp, query =="ZK218.8" &length > 100)
  if (nrow(df_left_8)>0 & nrow(df_1st_gap)>0){
    check_dir <- df_1st_gap[which.max(df_left_8$length),]
    if (check_dir$sstart > check_dir$send){
      left_df[nrow(left_df)+1,] <- c(readName,(min(df_1st_gap$send)- min(df_5s_1$send)), "-")
    }else{
      left_df[nrow(left_df)+1,] <- c(readName,(max(df_5s_1$send)-max(df_1st_gap$send)),"+")
    }
  } else if (nrow(df_1st_gap)>0){
    check_dir <- df_1st_gap[which.max(df_1st_gap$length),]
    if (check_dir$sstart > check_dir$send){
      left_df[nrow(left_df)+1,] <- c(readName,(check_dir$send- min(df_5s_1$send)), "-")
    }else{
      left_df[nrow(left_df)+1,] <- c(readName,(max(df_5s_1$send)-check_dir$send), "+")
    }
  }
}
#remove unspecific reads
left_df$remain =as.numeric(left_df$remain)
left_df_filter <- subset(left_df, remain >1000)
#sort
plot_length <- left_df_filter[with(left_df_filter, order(strand,-remain)),]
F_reads <- subset(plot_length, strand == "+")
R_reads <- subset(plot_length, strand == "-")
F_reads$idx <- c(seq(1:nrow(F_reads)))
R_reads$idx <- c(seq(1:nrow(R_reads)))
p<-ggplot()+
  scale_x_continuous(expand = c(0, 0),breaks=seq(0,67000,by=5000),labels=scales::comma,limits = c(0,67000))+
  scale_y_continuous(expand = c(0, 0),limits = c(-240,240),breaks=seq(-240,240,by=10))+
  geom_rect(data=F_reads,mapping=aes(xmin=0,xmax=remain,ymin=idx,ymax = idx+1),fill = "black")+
  geom_rect(data=R_reads,mapping=aes(xmin=0,xmax=remain,ymin=-idx,ymax = -idx-1),fill = "black")+
  geom_hline(yintercept=0,color= "red")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
p
pdf("5S_flanking_depth/N2YA_left.pdf",width=10,height = 8)
print(p)
dev.off()



##################N2YA_5s_only######################
rdna_list <- read.table("plot/N2YA/5s_only.list")
rdna_df <- data.frame(readID = character(), 
                      qstart =integer(),
                      qend= integer(),
                      strand = character(), 
                      stringsAsFactors = FALSE)
for (rdna_read in rdna_list$V1){
  readName <-as.character(rdna_read)
  df_tmp <-subset(df,df$subject == readName)
  df_5s_1 <-subset(df_tmp, query =="cel-5s-unit1" & length >42)
  check_dir <-df_5s_1 [which.max(df_5s_1$length),]
  if (check_dir$sstart > check_dir$send){
    align_length <-max(df_5s_1$sstart)-min(df_5s_1$send)
    rdna_df[nrow(rdna_df)+1,] <- c(readName,-align_length/2, align_length/2,"-")
  }else{
    align_length <-max(df_5s_1$send)-min(df_5s_1$sstart)
    rdna_df[nrow(rdna_df)+1,] <- c(readName,-align_length/2, align_length/2,"+")
  }
}
rdna_df$qstart=as.numeric(rdna_df$qstart)
rdna_df$qend=as.numeric(rdna_df$qend)
#sort
F_reads <- subset(rdna_df, strand == "+" & qend >1000)
R_reads <- subset(rdna_df, strand == "-" & qend >1000)
F_reads$sum2 <- F_reads$qend-F_reads$qstart
R_reads$sum2 <- R_reads$qend-R_reads$qstart

F_reads_sort <- F_reads[with(F_reads,order(-sum2)),]
F_reads_sort$idx <-c(seq(1:nrow(F_reads_sort)))

R_reads_sort <- R_reads[with(R_reads,order(-sum2)),]
R_reads_sort$idx <-c(seq(1:nrow(R_reads_sort)))

p<-ggplot()+
  scale_x_continuous(expand = c(0, 0),
                     breaks=seq(-33500,33500,by=5000),
                     labels=scales::comma,
                     limits = c(-33500,33500))+
  scale_y_continuous(expand = c(0, 0),limits = c(-240,240),breaks=seq(-240,240,by=10))+
  geom_rect(data=F_reads_sort,mapping=aes(xmin=qstart,xmax=qend,ymin=idx,ymax = idx+1),fill = "black")+
  geom_rect(data=R_reads_sort,mapping=aes(xmin=qstart,xmax=qend,ymin=-idx,ymax = -idx-1),fill = "black")+
  geom_hline(yintercept=0,color= "red")+
  geom_vline(xintercept=0,color= "red")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
p
pdf("5S_flanking_depth/N2YA_5s_only.pdf",width=10,height = 8)
print(p)
dev.off()

##################N2YA_right######################
right_list <- read.table("plot/N2YA/right_read.list")
right_df <- data.frame(readID = character(), 
                       remain= integer(),
                       strand = character(), 
                       stringsAsFactors = FALSE)
for (right_read in right_list$V1){
  readName <-as.character(right_read)
  df_tmp <-subset(df,df$subject == readName)
  df_5s_1 <-subset(df_tmp, query =="cel-5s-unit1" & length >42)
  df_variant <- subset(df_tmp, query =="corrected_5S_unit3:750-815" & qstart <20 & qend >45)
  
  if (nrow(df_variant)>0){
    check_dir <- df_variant[which.max(df_variant$length),]
    if (check_dir$sstart > check_dir$send){
      right_df[nrow(right_df)+1,] <- c(readName,-max(df_5s_1$sstart)+max(df_variant$sstart), "-")
    }else{
      right_df[nrow(right_df)+1,] <- c(readName,-min(df_variant$sstart),"+")
    }
  } 
  
}
#remove unspecific reads
right_df$remain = as.numeric(right_df$remain)
right_df_filter <- subset(right_df, remain <= -1000)

right_df_filter$remain=as.numeric(right_df_filter$remain)
#sort
plot_length <- right_df_filter[with(right_df_filter, order(strand,remain)),]

F_reads <- subset(plot_length, strand == "+")
R_reads <- subset(plot_length, strand == "-")
F_reads$idx <- c(seq(1:nrow(F_reads)))
R_reads$idx <- c(seq(1:nrow(R_reads)))
p<-ggplot()+
  scale_x_continuous(expand = c(0, 0),breaks=seq(-67000,0,by=5000),labels=scales::comma,limits = c(-67000,0))+
  scale_y_continuous(expand = c(0, 0),limits = c(-240,240),breaks=seq(-240,240,by=10))+
  geom_rect(data=F_reads,mapping=aes(xmin=0,xmax=remain,ymin=idx,ymax = idx+1),fill = "black")+
  geom_rect(data=R_reads,mapping=aes(xmin=0,xmax=remain,ymin=-idx,ymax = -idx-1),fill = "black")+
  geom_hline(yintercept=0,color= "red")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
p
pdf("5S_flanking_depth/N2YA_right.pdf",width=10,height = 8)
print(p)
dev.off()


##################N2YA_5s_only_readlength_distribution######################
rdna_list <- read.table("plot/N2YA/5s_only.list")
rdna_df <- data.frame(readID = character(), 
                      qstart =integer(),
                      qend= integer(),
                      strand = character(), 
                      stringsAsFactors = FALSE)
for (rdna_read in rdna_list$V1){
  readName <-as.character(rdna_read)
  df_tmp <-subset(df,df$subject == readName)
  df_5s_1 <-subset(df_tmp, query =="cel-5s-unit1" & length >42)
  check_dir <-df_5s_1 [which.max(df_5s_1$length),]
  if (check_dir$sstart > check_dir$send){
    align_length <-max(df_5s_1$sstart)-min(df_5s_1$send)
    rdna_df[nrow(rdna_df)+1,] <- c(readName,0, align_length,"-")
  }else{
    align_length <-max(df_5s_1$send)-min(df_5s_1$sstart)
    rdna_df[nrow(rdna_df)+1,] <- c(readName,0, align_length,"+")
  }
}
rdna_df$qstart=as.numeric(rdna_df$qstart)
rdna_df$qend=as.numeric(rdna_df$qend)
#sort
F_reads <- subset(rdna_df, strand == "+" & qend >1000)
R_reads <- subset(rdna_df, strand == "-" & qend >1000)
F_reads$sum2 <- F_reads$qend-F_reads$qstart
R_reads$sum2 <- R_reads$qend-R_reads$qstart

F_reads_sort <- F_reads[with(F_reads,order(-qend)),]
F_reads_sort$idx <-c(seq(1:nrow(F_reads_sort)))

R_reads_sort <- R_reads[with(R_reads,order(-qend)),]
R_reads_sort$idx <-c(seq(1:nrow(R_reads_sort)))

p<-ggplot()+
  scale_x_continuous(expand = c(0, 0),
                     breaks=seq(0,67000,by=5000),
                     labels=scales::comma,
                     limits = c(0,67000))+
  scale_y_continuous(expand = c(0, 0),limits = c(-240,240),breaks=seq(-240,240,by=10))+
  geom_rect(data=F_reads_sort,mapping=aes(xmin=qstart,xmax=qend,ymin=idx,ymax = idx+1),fill = "black")+
  geom_rect(data=R_reads_sort,mapping=aes(xmin=qstart,xmax=qend,ymin=-idx,ymax = -idx-1),fill = "black")+
  geom_hline(yintercept=0,color= "red")+
  geom_vline(xintercept=0,color= "red")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
p
pdf("5S_flanking_depth/N2YA_5s_only_readlength.pdf",width=10,height = 8)
print(p)
dev.off()
