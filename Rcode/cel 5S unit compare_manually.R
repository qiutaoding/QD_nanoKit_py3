library("ggplot2")
p<-ggplot()+
  xlab("cel-5s unit2: 946 bp")+
  ylab("cel-5S unit1: 976 bp")+
  scale_x_continuous(breaks=seq(0,1000,by=500),limits=c(0,1000),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,1000,by=500),limits=c(0,1000),expand=c(0,0))+
  geom_segment(aes(x=946,y=0,xend=946,yend=976),color="red",linetype="dashed",size=0.5)+
  geom_segment(aes(x=0,y=976,xend=946,yend=976),color="red",linetype="dashed",size=0.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=10))+
  geom_segment(mapping = aes(x=0,y=0,xend=791,yend=791),size=1,color= "black",linetype=1)+
  geom_segment(mapping = aes(x=775,y=805,xend=946,yend=976),size=1,color= "black",linetype=1)
p
setwd("D:/16483960_drive/Project/Project2_rDNA/paper figure")
pdf("cel-5S unit1-cel-5S unit 2 dotplot.pdf",width=6,height = 6)
print(p)
dev.off()

