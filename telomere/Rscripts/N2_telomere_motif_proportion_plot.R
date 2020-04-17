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
binSet = 1
#start plot
#The motif coverage and proportion in each strand of chromosome end has two plots in one PDF
for (i in stage){
  for (j in chr){
    for (k in position) {
      for (l in strand) {
        Cov_df = read.table(paste0( "tidehunter/stats/coverage/", i , "_", j , k , "_" , l , "_Coverage.tsv"),  header = TRUE)
        motif_df = read.table(paste0("tidehunter/stats/sorted/" , i , "_", j , k , "_" , l , "_stats_sorted.tsv"),  header = TRUE , stringsAsFactors = FALSE)
        #only select the motif with >= 200 copies
        motif_df <- subset(motif_df , CopyNum >= 200)
        #The canonical telomere repeat motif TTAGGC is not always been the maximum copy number
        #use which() function to locate the number row iwth canonical repeat
        motif_df[which(motif_df[,1] =="canonical"),1] <- "TTAGGC" #impportant to convert RepeatSeq name, otherwise error occurs
        Cov_df <-subset(Cov_df , select = motif_df[,1])#subset df by the RepeatSeq in motif_df name 
        for (x  in 1:ncol(Cov_df)){
          aa <- subset(motif_df, RepeatSeq == colnames(Cov_df)[x])$Consensus#Get motif sequences from motif df
          #colnames(Cov_df)[x] <- aa
          #convert telomere repeat motif consistent
          #assign motif sequences
          if (aa %in% TTGGGC){
            colnames(Cov_df)[x] <- "TTGGGC"
          }else if ( aa %in% TTAGC){
            colnames(Cov_df)[x] <- "TTAGC"
          }else if ( aa %in% TCAGGC){
            colnames(Cov_df)[x] <- "TCAGGC"
          }else if ( aa %in% TAGGC){
            colnames(Cov_df)[x] <- "TAGGC"
          }else if ( aa %in% AGGCC){
            colnames(Cov_df)[x] <- "AGGCC"
          }else if ( aa %in% AGGC){
            colnames(Cov_df)[x] <- "AGGC"
          }else{colnames(Cov_df)[x] <- aa
          }
        }
        #add location information to df, default is 10 kb df
        Cov_df <- cbind(position = seq(1 , 10000) , Cov_df)
        #melt repeat df suitable for plotting
        mt_Cov_df <- melt(Cov_df , id.vars = c("position"), variable.name = "Repeat", value.name = "Depth")
        tmp_plot1 <- ggplot(mt_Cov_df, aes(x = position, y = Depth, fill = Repeat )) + 
          geom_bar(position="fill", stat="identity")+
          scale_x_continuous(breaks=seq(0 , 10000 , by = 500),
                             expand = c(0 , 0))+  
          scale_y_continuous(expand = c(0 , 0),
                             breaks=seq(0 , 1 ,by = 0.1),#y-axis breaks
                             labels= seq(0 , 100,by = 10))+ 
          coord_cartesian(xlim = c(0 , 5500))+
          ylab("Proportion %")+
          theme(panel.grid.major = element_blank(),
                axis.title.x=element_blank(),#remove x label
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.text = element_text(size=10, colour = "black"),
                axis.title = element_text(size=12, colour = "black"))
        tmp_plot2 <- ggplot(mt_Cov_df ,aes( x = position, weight = Depth / binSet , fill = Repeat ))+
          geom_histogram(aes( x = position, weight = Depth / binSet , fill = Repeat ),binwidth = binSet )+
          scale_x_continuous(breaks=seq(0 , 10000 , by = 500),
                             expand = c(0 , 0))+  
          scale_y_continuous(expand = c(0 , 0),
                             labels= seq(0 , 10 * (as.integer(max(mt_Cov_df$Depth)/10) + 1) ,by = 10),
                             breaks=seq(0 ,  10 * (as.integer(max(mt_Cov_df$Depth)/10) + 1) ,by = 10))+#y-axis breaks 
          coord_cartesian(xlim = c(0 , 5500))+
          ylab("Coverage")+
          theme(panel.grid.major = element_blank(),
                axis.title.x=element_blank(),#remove x label
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.text = element_text(size=10, colour = "black"),
                axis.title = element_text(size=12, colour = "black"))
        
        pdf(paste0("plot_file/motif_proportion/P8p4k3/", i , "_", j , k , "_" , l ,"_proportion_6X6.pdf") , width = 6 , height = 6)
        grid.arrange( tmp_plot2, tmp_plot1 , nrow = 2)
        dev.off()
        print(paste("====",i , j , k ,l , "finished===="))
      }
    }
  }
}
rm(i,j,k,l, tmp_plot1, tmp_plot2, aa)
