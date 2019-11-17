setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/plot/")
library("ggplot2")
df_N2EM=read.table("N2EM_5s_only_read.fa.fai",header=FALSE)
df_N2L1=read.table("N2L1_5s_only_read.fa.fai",header=FALSE)
df_N2YA=read.table("N2YA_5s_only_read.fa.fai",header=FALSE)
df_zzy0600=read.table("zzy0600_5s_only_read.fa.fai",header=FALSE)
df_zzy0603=read.table("zzy0603_5s_only_read.fa.fai",header=FALSE)

p_des<-ggplot()+
  geom_density(data=df_N2EM,aes(x=V2),color = "red")+
  geom_density(data=df_N2L1,aes(x=V2),color = "blue")+
  geom_density(data=df_N2YA,aes(x=V2),color = "green")+
  geom_density(data=df_zzy0600,aes(x=V2),color = "pink")+
  geom_density(data=df_zzy0603,aes(x=V2),color = "violet")
p_des

library(ggpubr)
binset= 1000
p1<-ggplot()+geom_histogram(data=df_N2EM,aes(x=V2),color = "red", binwidth =binset )+
  scale_x_continuous(expand = c(0, 0),breaks=seq(10000,70000,by=10000),labels=scales::comma,limits = c(10000,70000))
p2<-ggplot()+geom_histogram(data=df_N2L1,aes(x=V2),color = "blue", binwidth =binset)+
  scale_x_continuous(expand = c(0, 0),breaks=seq(10000,70000,by=10000),labels=scales::comma,limits = c(10000,70000))
p3<-ggplot()+geom_histogram(data=df_N2YA,aes(x=V2),color = "green", binwidth =binset)+
  scale_x_continuous(expand = c(0, 0),breaks=seq(10000,70000,by=10000),labels=scales::comma,limits = c(10000,70000))
p4<-ggplot()+geom_histogram(data=df_zzy0600,aes(x=V2),color = "pink", binwidth =binset)+
  scale_x_continuous(expand = c(0, 0),breaks=seq(10000,70000,by=10000),labels=scales::comma,limits = c(10000,70000))
p5<-ggplot()+geom_histogram(data=df_zzy0603,aes(x=V2),color = "violet", binwidth =binset)+
  scale_x_continuous(expand = c(0, 0),breaks=seq(10000,70000,by=10000),labels=scales::comma,limits = c(10000,70000))
ggarrange(p1, p2,p3,p4,p5, nrow = 5)

p5
