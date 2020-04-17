setwd("F:/nanopore_data/telomere/N2_extract/trimming_result")
library("ggplot2")
####load telomere sequences length and merge to single df####
stage <- c("N2EM","N2L1","N2YA")
chr <-  c("I", "II","III","IV","V","X")
position <- c("L","R")
strand <- c("F" , "R")
fai_title <- c("readname", "length","offset","linebase","linewidth")
#strand specific df
df_ONT_sum = data.frame(stage = character(),#library from C.elegans developmental stage, e.g. N2EM
                        chromosome = character(),#chromosome, e.g. I
                        position = character(),#chromosome end, e.g.  L
                        strand = character(), #sequence strand, e.g. F
                        readname = character(),#nanopore reads read name
                        length = integer(),#telomere sequence length
                        stringsAsFactors = FALSE)
#pre-request bash command: for f in trimming_result/strand/*.fa;do samtools faidx $f; done
for (i in stage){
  for (j in chr){
    for (k in position) {
      for (l in strand){
        input_df = read.table(paste0( "ONT_reads/" , i , "_" , j , k , "_" , l , ".fa.fai"), header = FALSE)#read fasta index file from samtools
        colnames(input_df) <- fai_title#add column name
        input_df <- cbind(stage = i , chromosome = j, position = k , strand = l, input_df)#add stage, chromosome, position, strand information before index df
        input_df <-subset(input_df , select = -c(offset,linebase , linewidth))#remove unnecessary columns from df
        df_ONT_sum <- rbind(df_ONT_sum , input_df)#append read length information to the main df for further analysis and plotting
      }
    }
  }
}
rm(input_df , i, j ,k , l)
#write merged df to file to save processed data
write.csv(df_ONT_sum,file = "plot_file/ONT_reads_with_telomere_length_strand_table.csv", row.names=FALSE)#save file, sep by comma, and keep header

#df_ONT_sum = read.table("plot_file/ONT_reads_with_telomere_length_strand_table.csv" , header = TRUE , sep =",")

max_label = 9000#longest read length
max_zoom = 5000#maximum of visual zoom
break_by = 1000#read length breaks

#######sep by stage only########
ca <- ggplot(df_ONT_sum, aes(x = stage , y = length , fill = stage))+
  geom_violin()+#add normal violin plot layer
  geom_boxplot(width=0.2 , notch = TRUE , 
               outlier.shape = NA,
               position = position_dodge()
               )+#add box plot, add 95% confidence notch, set outliner size
  scale_y_continuous(breaks=seq(0 , 100000 ,by = 20000),#y-axis breaks
                     labels= seq(0 ,100,by = 20),#rename y-axis break labels
                     expand = c(0 , 0))+#set plot to 0,0 poistion
  coord_cartesian(ylim = c(0 , 100000))+#show visual zoom 
  #facet_grid(chromosome ~ position)+#arrange plot by chromosome and chromosome end
  #set legend format
  scale_fill_discrete(name="Stage",
                      breaks=stage,#c("N2EM", "N2L1", "N2YA")
                      labels=c("EMB", "L1", "YA"))+
  #tick mark labels
  scale_x_discrete(labels=c("N2EM" = "EMB",
                            "N2L1" = "L1",
                            "N2YA" = "YA"))+
  #X-axis label
  #xlab(expression(paste(italic("C. elegans")," N2 chromosome ends")))+
  ylab("Nanopore reads length (kb)")+
  theme(axis.title.x=element_blank(),#remove x label
        panel.grid.major = element_blank(),#remove major grid
        panel.grid.minor = element_blank(),#remove minor grid
        panel.background = element_blank(), #remove background color
        axis.line = element_line(colour = "black"),#set axis-line colour
        legend.position = "none",#remove legend
        #legend.position="top",
        #legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),#set axis text size and colour
        axis.title = element_text(size=12, colour = "black"))#set axis title size and colour
ca
#export as pdf 3X5 inches
aov_stage<- aov(df_ONT_sum$length ~ df_ONT_sum$stage)
summary(aov_stage)
tuk_stage <- TukeyHSD(aov_stage )
tuk_stage



#############sep by strand#############
cb <- ggplot(df_ONT_sum, aes(x = strand,y = length, fill = strand))+
  geom_violin()+#add normal violin plot layer
  geom_boxplot(width=0.2 , notch = TRUE , 
               outlier.shape = NA,
               position = position_dodge()
  )+#add box plot, add 95% confidence notch, set outliner size
  scale_y_continuous(breaks=seq(0 , 100000 ,by = 20000),#y-axis breaks
                     labels= seq(0 ,100,by = 20),#rename y-axis break labels
                     expand = c(0 , 0))+#set plot to 0,0 poistion
  coord_cartesian(ylim = c(0 , 100000))+#show visual zoom 
  #facet_grid(chromosome ~ position)+#arrange plot by chromosome and chromosome end
  #set legend format
  scale_fill_discrete(name="Stage",
                      breaks=stage,#c("N2EM", "N2L1", "N2YA")
                      labels=c("EMB", "L1", "YA"))+
  #tick mark labels
  scale_x_discrete(labels=c("N2EM" = "EMB",
                            "N2L1" = "L1",
                            "N2YA" = "YA" ,
                            "F" = "G-strand",
                            "R" = "C-strand"))+
  #X-axis label
  #xlab(expression(paste(italic("C. elegans")," N2 chromosome ends")))+
  ylab("Nanopore reads length (kb)")+
  theme(axis.title.x=element_blank(),#remove x label
        panel.grid.major = element_blank(),#remove major grid
        panel.grid.minor = element_blank(),#remove minor grid
        panel.background = element_blank(), #remove background color
        axis.line = element_line(colour = "black"),#set axis-line colour
        legend.position = "none",#remove legend
        #legend.position="top",
        #legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),#set axis text size and colour
        axis.title = element_text(size=12, colour = "black"))#set axis title size and colour
cb
#export as pdf 3X5 inches
aov_strand<- aov(df_ONT_sum$length ~ df_ONT_sum$strand)
summary(aov_strand)
tuk_strand <- TukeyHSD(aov_strand )
tuk_strand
