setwd("F:/nanopore_data/telomere/N2_extract/trimming_result")
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
stage <- c("N2EM","N2L1","N2YA")
chr <-  c("I", "II","III","IV","V","X")
position <- c("L","R")
strand <- c("F" , "R")
enumerateRepeat <- function(motif){
  motif_element <- strsplit(motif, "")#split characters in telomere repeat motif
  motif_vector <- c()#make an empty vector for return
  for (i in 1:length(motif_element[[1]])){
    tmp1 <- paste0(motif_element[[1]][i:length(motif_element[[1]])] , collapse ="")
    tmp2 <- paste0(motif_element[[1]][0:(i-1)] , collapse ="")
    motif_vector <- c(motif_vector, paste0(tmp1 , tmp2 ))
  }
  return(motif_vector)
}
TTAGGC <- enumerateRepeat("TTAGGC")
TTGGGC <- enumerateRepeat("TTGGGC")
TTAGC <- enumerateRepeat("TTAGC")
TCAGGC <- enumerateRepeat("TCAGGC")
TAGGC <- enumerateRepeat("TTAGC")
AGGCC <- enumerateRepeat("AGGCC")
AGGC <- enumerateRepeat("AGGC")
tide_title <- c("readName", "consN","readLen","start","end" ,
                "consLen", "copyNum", "fullLen", "subPos", "consensus")
colorPalette <- c("TTAGGC" = "#2c7bb6",
                  "TTGGGC" = "#fdae61",
                  "AGGC" = "#ffffbf",
                  "Others" = "#abd9e9")
for (i in stage){
  for (j in chr){
    for (k in position){
      for (l in strand){
        tideInfo = read.table(paste0( "tidehunter/P8p4k3/", i , "_", j , k , "_" , l , "_P8p4k3.tsv"),  header = FALSE, stringsAsFactors = FALSE)
        colnames(tideInfo) <- tide_title
        #df sorted by descending telomere sequence length
        order_tideInfo <- tideInfo[order(-tideInfo$readLen), ]
        rm(tideInfo)
        #convert telomere motif repeat in the ordered df
        for ( x in 1:nrow(order_tideInfo)){
          aa <- order_tideInfo[x,]$consensus
          if(aa %in% TTAGGC){
            order_tideInfo[x,]$consensus <- "TTAGGC"
          }else if (aa %in% TTGGGC){
            order_tideInfo[x,]$consensus <- "TTGGGC"
          }else if (aa %in% AGGC){
            order_tideInfo[x,]$consensus <- "AGGC"
          }else{order_tideInfo[x,]$consensus <- "Others"}
        }
        uniq_list <- unique(order_tideInfo$readName)
        for ( x in 1:nrow(order_tideInfo)){
          order_tideInfo[x,]$readName <- which(uniq_list == order_tideInfo[x,]$readName)
        }
        order_tideInfo$readName = as.numeric(order_tideInfo$readName)
        #initialize ggplot 
        pp <- ggplot(data = order_tideInfo) + 
          geom_rect(aes(xmin = start, xmax = end, 
                        ymin =  readName - 1, 
                        ymax = readName,
                        fill = consensus))+
          scale_fill_manual(values = colorPalette,
                            name = "Motif")+#legend correction
          scale_x_continuous(breaks=seq(0 , 10000 , by = 500),
                             expand = c(0 , 0))+  
          scale_y_continuous(expand = c(0 , 0),
                             labels = seq(0 , 200 ,by = 10),
                             breaks = seq(0 , 200 ,by = 10))+
          coord_cartesian(xlim = c(0 , 5500))+
          ylab("Telomere repeat motifs in each read")+
          theme(panel.grid.major = element_blank(),
                axis.title.x=element_blank(),#remove x label
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.text = element_text(size=10, colour = "black"),
                axis.title = element_text(size=12, colour = "black"))
        pdf(paste0("plot_file/motif_cluster_read/", i , "_", j , k , "_" , l ,"_P8p4k3_6X6.pdf") ,
            width = 6 , 
            height = 6)
        print(pp)
        dev.off()
        print(paste("====",i , j , k ,l , "finished===="))
      }
    }
  }
}

rm(aa, pp, i,j,k,l,x)
