
library("ggplot2")
library("ggtree")
setwd("D:/16483960_drive/Project/Project2_rDNA/version/V5")

ratio <- read.csv("5Sunit_330strain.txt",sep="\t",header =T)
colnames(ratio)<-c("strain","unit1","unit2","percent")

tree <- read.tree("330_genome.raxml.bestTree")
p1<- ggtree(tree)
p1
#add label
p2<- p1+geom_tiplab()
p2
#rotate
p3 <-open_tree(p1,180)
p3

#circular without branch length
p4<- ggtree(tree,
            layout="circular",
            branch.length = "none")
 
p4
#circular looks fine
p5<- p4 + geom_tiplab(size= 1,aes(angle=angle))
p5
#add color test >1%
groupInfo <- ratio[,1][1:157]
with_unit2 <-groupOTU(tree,groupInfo)
p6<- ggtree(with_unit2, 
            aes(color=group),
            layout="circular",
            branch.length = "none")+ 
  #scale_color_manual(values=c("pink", "blue"))+
  geom_tiplab(size= 1,aes(angle=angle))
p6
#add color test >1% and manually select strain >0.1% and reads no. >10
groupInfo <- ratio[,1][1:157]
groupInfo <-as.character(groupInfo)
groupInfo <-c(groupInfo,'JU1568',"MY2535","JU1580","ECA778","JU782","CX11264","NIC513")
with_unit2 <-groupOTU(tree,groupInfo)
p6<- ggtree(with_unit2, 
            aes(color=group),
            layout="circular",
            branch.length = "none")+ 
  #scale_color_manual(values=c("pink", "blue"))+
  geom_tiplab(size= 1,aes(angle=angle))
p6




