#set working directry
setwd("F:/nanopore_data/rDNA/recall_rdna_analysis/SNP/for_new_variant/")
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
#preset data 
cel <-c("N2EM", "N2L1", "N2YA", "CB4856")
cbr <- c("AF16")
cel_rdna <- c("5S")
cbr_rdna <-c("5S1", "5S2")
coverage = c(N2EM= 100.272, N2L1 = 35.6229, N2YA = 48.2667, AF16=90.8984, CB4856 = 128.261)
seg_hight = 9
txt_hight = 9.5

#modify depends on demand, plot_indel() function should be run before run the plotting command
#5S rDNA
plot_indel( "N2EM", "5S" )
plot_indel( "N2L1", "5S" )
plot_indel( "N2YA", "5S" )
plot_indel( "CB4856", "5S" )
plot_indel( "AF16", "5S" )
plot_indel( "AF16", "5S2" )
#45S rDNA
plot_45S_indel( "N2EM", "45S" )
plot_45S_indel( "N2L1", "45S" )
plot_45S_indel( "N2YA", "45S" )
plot_45S_indel( "CB4856", "45S" )
plot_45S_indel( "AF16", "45S" )

#tables were generated from SAM files with read_SAM_for_INDEL.ipynb 
plot_indel <- function(stage , rDNA){
  #read in processed data
  df2=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 2 ,".tsv"), header =T , sep = "\t")
  df5=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 5 ,".tsv"), header =T , sep = "\t")
  df10=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 10 ,".tsv"), header =T , sep = "\t")
  #rRNA and sl1 start and end position
  if (stage %in% cel){
    rRNA_start <- 177
    rRNA_end <- 295
    SL1_start <-976
    SL1_end <-879
  }else if (stage %in% cbr){
    if (rDNA =="5S"){
      rRNA_start <- 218
      rRNA_end <- 336
      SL1_start <-938
      SL1_end <- 842
    }else if (rDNA == "5S2"){
      rRNA_start <- 171
      rRNA_end <- 289
      SL1_start <-600
      SL1_end <- 696
    }
  }
  #INDEL size > 10 in ONT reads mapping results
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
  #INDEL size > 5 in ONT reads mapping results
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
  #INDEL size > 2 in ONT reads mapping results
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
  #print three plot in single PDF
  pdf(paste0("plot/", stage,"_", rDNA, "_bins.pdf"), width = 8, height = 12) # Open a new pdf file
  grid.arrange(p2, p5,p10, nrow = 3) # Write the grid.arrange in the file
  dev.off()
}

plot_45S_indel <- function(stage , rDNA){
  df2=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 2 ,".tsv"), header =T , sep = "\t")
  df5=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 5 ,".tsv"), header =T , sep = "\t")
  df10=read.table(paste0("indel_df/",stage,"_", rDNA,"_bin_", 10 ,".tsv"), header =T , sep = "\t")
  if(stage %in% cel){      
    rRNA_start.18S <- 934
    rRNA_end.18S <- 2687
    rRNA_start.5_8S <- 3152
    rRNA_end.5_8S <- 3304
    rRNA_start.26S <- 3689
    rRNA_end.26S <- 7197
  }else if(stage %in% cbr){
    rRNA_start.18S <- 1263
    rRNA_end.18S <- 3016
    rRNA_start.5_8S <- 3460
    rRNA_end.5_8S <- 3612
    rRNA_start.26S <- 3984
    rRNA_end.26S <- 7487
  }
  p10 <- ggplot(df10) + 
    geom_histogram( aes(pos, weight = DEL/coverage[stage]),
                    binwidth = 1 , fill = "#F8766D")+
    geom_histogram( aes(pos, weight = INS/coverage[stage]),
                    binwidth = 1, fill = "#619CFF")+
    coord_cartesian(xlim= c(0,7600),
                    ylim=c(0,10))+
    scale_x_continuous(label =comma, breaks = seq(0, 10000, by = 1000),
                       expand=c(0,0))+
    scale_y_continuous(label =seq(0, 20, by = 2), breaks = seq(0, 20, by = 2),
                       expand=c(0,0))+
    ylab("Copy number: INDEL>10")+
    geom_segment(aes(x = rRNA_start.18S , y= seg_hight , xend = rRNA_end.18S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.18S , rRNA_end.18S)) ,
             y= txt_hight, label = "18S rRNA", color = "black")+
    geom_segment(aes(x = rRNA_start.5_8S , y= seg_hight , xend = rRNA_end.5_8S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.5_8S , rRNA_end.5_8S)) ,
             y= txt_hight, label = "5.8S rRNA", color = "black")+
    geom_segment(aes(x = rRNA_start.26S , y= seg_hight , xend = rRNA_end.26S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.26S , rRNA_end.26S)) ,
             y= txt_hight, label = "26S", color = "black")+
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
    coord_cartesian(xlim= c(0,7600),
                    ylim=c(0,10))+
    scale_x_continuous(label =comma, breaks = seq(0, 10000, by = 1000),
                       expand=c(0,0))+
    scale_y_continuous(label =seq(0, 20, by = 2), breaks = seq(0, 20, by = 2),
                       expand=c(0,0))+
    ylab("Copy number: INDEL>10")+
    geom_segment(aes(x = rRNA_start.18S , y= seg_hight , xend = rRNA_end.18S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.18S , rRNA_end.18S)) ,
             y= txt_hight, label = "18S rRNA", color = "black")+
    geom_segment(aes(x = rRNA_start.5_8S , y= seg_hight , xend = rRNA_end.5_8S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.5_8S , rRNA_end.5_8S)) ,
             y= txt_hight, label = "5.8S rRNA", color = "black")+
    geom_segment(aes(x = rRNA_start.26S , y= seg_hight , xend = rRNA_end.26S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.26S , rRNA_end.26S)) ,
             y= txt_hight, label = "26S", color = "black")+
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
    coord_cartesian(xlim= c(0,7600),
                    ylim=c(0,10))+
    scale_x_continuous(label =comma, breaks = seq(0, 10000, by = 1000),
                       expand=c(0,0))+
    scale_y_continuous(label =seq(0, 20, by = 2), breaks = seq(0, 20, by = 2),
                       expand=c(0,0))+
    ylab("Copy number: INDEL>10")+
    geom_segment(aes(x = rRNA_start.18S , y= seg_hight , xend = rRNA_end.18S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.18S , rRNA_end.18S)) ,
             y= txt_hight, label = "18S rRNA", color = "black")+
    geom_segment(aes(x = rRNA_start.5_8S , y= seg_hight , xend = rRNA_end.5_8S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.5_8S , rRNA_end.5_8S)) ,
             y= txt_hight, label = "5.8S rRNA", color = "black")+
    geom_segment(aes(x = rRNA_start.26S , y= seg_hight , xend = rRNA_end.26S , yend = seg_hight))+
    annotate(geom="text" , x= mean(c(rRNA_start.26S , rRNA_end.26S)) ,
             y= txt_hight, label = "26S", color = "black")+
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


