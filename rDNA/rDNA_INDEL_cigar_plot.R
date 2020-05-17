setwd("D:/rDNA/")
library("ggplot2")
library(scales)
library(grid)
library(gridExtra)
#preset
cel <-c("N2EM", "N2L1", "N2YA", "CB4856")
cbr <- c("AF16")
cel_rdna <- c("5S")
cbr_rdna <-c("5S_unit1", "5S_unit2")
coverage = c(N2EM= 100.272, N2L1 = 35.6229, N2YA = 48.2667, AF16=90.8984, CB4856 = 128.261)
seg_hight = 9
txt_hight = 9.5

#modify depends on demand
stage =  "N2L1"
rDNA =  "5S"
plot_indel( "AF16", "5S2" )

#tables were generated from SAM files with read_SAM_for_INDEL.ipynb 
plot_indel <- function(stage , rDNA){
  df2=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 2 ,".tsv"), header =T , sep = "\t")
  df5=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 5 ,".tsv"), header =T , sep = "\t")
  df10=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 10 ,".tsv"), header =T , sep = "\t")
  
  if (stage %in% cel){
    rRNA_start <- 177
    rRNA_end <- 295
    SL1_start <-976
    SL1_end <-879
  }else if (stage %in% cbr_rdna){
    if (rDNA =="5S_unit1"){
      rRNA_start <- 218
      rRNA_end <- 336
      SL1_start <-938
      SL1_end <- 842
    }else if (rDNA == "5S_unit2"){
      rRNA_start <- 169
      rRNA_end <- 287
      SL1_start <-597
      SL1_end <- 693
    }
  }
  p10 <- ggplot(df10) + 
    geom_histogram( aes(pos, weight = DEL/coverage[stage]),
                    binwidth = 1 , fill = "#F8766D")+
    geom_histogram( aes(pos, weight = INS/coverage[stage]),
                    binwidth = 1, fill = "#619CFF")+
    coord_cartesian(xlim= c(0,1000),
                    ylim=c(0,10))+
    scale_x_continuous(label =comma, breaks = seq(0, 10000, by = 100),
                       expand=c(0,0))+
    scale_y_continuous(label =seq(0, 20, by = 2), breaks = seq(0, 20, by = 2),
                       expand=c(0,0))+
    ylab("Copy number: INDEL>10")+
    geom_segment(aes(x = rRNA_start , y= seg_hight , xend = rRNA_end , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start,rRNA_end)) ,
             y= txt_hight, label = "5S rRNA", color = "black")+
    geom_segment(aes(x = SL1_start , y= seg_hight , xend = SL1_end , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(SL1_start,SL1_end)) ,
             y= txt_hight, label = "SL1", color = "black")+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size= 12,colour = "black"),
          axis.title.x = element_blank())
  p5 <- ggplot(df5) + 
    geom_histogram( aes(pos, weight = DEL/coverage[stage]),
                    binwidth = 1 , fill = "#F8766D")+
    geom_histogram( aes(pos, weight = INS/coverage[stage]),
                    binwidth = 1, fill = "#619CFF")+
    coord_cartesian(xlim= c(0,1000),
                    ylim=c(0,10))+
    scale_x_continuous(label =comma, breaks = seq(0, 10000, by = 100),
                       expand=c(0,0))+
    scale_y_continuous(label =seq(0, 20, by = 2), breaks = seq(0, 20, by = 2),
                       expand=c(0,0))+
    ylab("Copy number: INDEL>5")+
    geom_segment(aes(x = rRNA_start , y= seg_hight , xend = rRNA_end , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start,rRNA_end)) ,
             y= txt_hight, label = "5S rRNA", color = "black")+
    geom_segment(aes(x = SL1_start , y= seg_hight , xend = SL1_end , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(SL1_start,SL1_end)) ,
             y= txt_hight, label = "SL1", color = "black")+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size= 12,colour = "black"),
          axis.title.x = element_blank())
  p2 <- ggplot(df2) + 
    geom_histogram( aes(pos, weight = DEL/coverage[stage]),
                    binwidth = 1 , fill = "#F8766D")+
    geom_histogram( aes(pos, weight = INS/coverage[stage]),
                    binwidth = 1, fill = "#619CFF")+
    coord_cartesian(xlim= c(0,1000),
                    ylim=c(0,10))+
    scale_x_continuous(label =comma, breaks = seq(0, 10000, by = 100),
                       expand=c(0,0))+
    scale_y_continuous(label =seq(0, 20, by = 2), breaks = seq(0, 20, by = 2),
                       expand=c(0,0))+
    ylab("Copy number: INDEL>2")+
    geom_segment(aes(x = rRNA_start , y= seg_hight , xend = rRNA_end , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start,rRNA_end)) ,
             y= txt_hight, label = "5S rRNA", color = "black")+
    geom_segment(aes(x = SL1_start , y= seg_hight , xend = SL1_end , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(SL1_start,SL1_end)) ,
             y= txt_hight, label = "SL1", color = "black")+
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size= 12,colour = "black"),
          axis.title.x = element_blank())
  
  
  pdf(paste0("plot/", stage,"_", rDNA, "_bins.pdf"), width = 8, height = 12) # Open a new pdf file
  grid.arrange(p2, p5,p10, nrow = 3) # Write the grid.arrange in the file
  dev.off()
}





