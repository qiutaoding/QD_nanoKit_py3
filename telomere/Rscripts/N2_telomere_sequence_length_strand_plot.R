#devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R")#source of flat violin
setwd("F:/nanopore_data/telomere/N2_extract/trimming_result")
library("ggplot2")
####load telomere sequences length and merge to single df####
stage <- c("N2EM","N2L1","N2YA")
chr <-  c("I", "II","III","IV","V","X")
position <- c("L","R")
strand <- c("F" , "R")
fai_title <- c("readname", "length","offset","linebase","linewidth")
#strand specific df
df_strand_sum = data.frame(stage = character(),#library from C.elegans developmental stage, e.g. N2EM
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
        input_df = read.table(paste0( "strand/" , i , "_" , j , k , "_" , l , ".fa.fai"), header = FALSE)#read fasta index file from samtools
        colnames(input_df) <- fai_title#add column name
        input_df <- cbind(stage = i , chromosome = j, position = k , strand = l, input_df)#add stage, chromosome, position, strand information before index df
        input_df <-subset(input_df , select = -c(offset,linebase , linewidth))#remove unnecessary columns from df
        df_strand_sum <- rbind(df_strand_sum , input_df)#append read length information to the main df for further analysis and plotting
      }
    }
  }
}
rm(input_df , i, j ,k , l)
#write merged df to file to save processed data
write.csv(df_strand_sum,file = "plot_file/telomere_length_strand_table.csv", row.names=FALSE)#save file, sep by comma, and keep header


df_strand_sum = read.table("plot_file/telomere_length_strand_table.csv" , header = TRUE , sep =",")

#set plot parameter
max_label = 9000#longest read length
max_zoom = 5000#maximum of visual zoom
break_by = 1000#read length breaks
#########start plotting###############
#########plot telomere length by strand ###############
a <- ggplot(df_strand_sum, aes(x = strand , y = length , fill = strand))+#by strand
  geom_violin(position = position_dodge(width = 1))+
  geom_boxplot(width = 0.2 , 
               notch = TRUE ,
               outlier.shape = NA, 
               position = position_dodge(width = 1))+
  #geom_jitter(width = 0.1 , alpha = 0.1)+
  scale_y_continuous(breaks=seq(0,max_label,by = break_by),#y-axis breaks
                     labels= seq(0,max_label/1000,by = break_by/1000),#rename y-axis break labels
                     expand = c(0, 0))+#set plot to 0,0 poistion
  coord_cartesian(ylim = c(0 , max_zoom))+#show visual zoom 
  #tick mark labels
  scale_x_discrete(labels = c("F" = "Forward","R" = "Reverse"))+
  #X-axis label
  #xlab(expression(paste(italic("C. elegans")," N2 chromosomes")))+
  ylab("Telomere sequences length (kb)")+
  #scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  #scale_fill_discrete(name = "Strand", labels = c("Forward", "Reverse"))+#legend correction
  theme(#panel.grid.major = element_blank(),
    axis.title.x=element_blank(),#remove x label
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.position = "none",#remove legend from plot
    #legend.position="top",
    #legend.direction = "horizontal",
    axis.text = element_text(size=10, colour = "black"),
    axis.title = element_text(size=12, colour = "black"))
a
#export as PDF 3X5 inches #1_N2_telomere_length_by_strand_3X5_box_violin.pdf

######### plot telomere length by strand and stage###############
b <- ggplot(df_strand_sum, aes(x = strand , y = length , fill = stage))+#by strand
  geom_violin(position = position_dodge(width = 1))+
  geom_boxplot(width = 0.2 , 
               notch = TRUE ,
               outlier.shape = NA, 
               position = position_dodge(width = 1))+
  #geom_jitter(width = 0.1 , alpha = 0.1)+
  scale_y_continuous(breaks=seq(0,max_label,by = break_by),#y-axis breaks
                     labels= seq(0,max_label/1000,by = break_by/1000),#rename y-axis break labels
                     expand = c(0, 0))+#set plot to 0,0 poistion
  coord_cartesian(ylim = c(0 , max_zoom))+#show visual zoom 
  #coord_cartesian(ylim = c(2000,2500))+#0 , max_zoom))+#show visual zoom 
  #tick mark labels
  scale_x_discrete(labels = c("F" = "Forward","R" = "Reverse"))+
  #X-axis label
  #xlab(expression(paste(italic("C. elegans")," N2 chromosomes")))+
  ylab("Telomere sequences length (kb)")+
  scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  #scale_fill_discrete(name = "Strand", labels = c("Forward", "Reverse"))+#legend correction
  theme(#panel.grid.major = element_blank(),
    axis.title.x=element_blank(),#remove x label
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    #legend.position = "none",#remove legend from plot
    #legend.position="top",
    #legend.direction = "horizontal",
    axis.text = element_text(size=10, colour = "black"),
    axis.title = element_text(size=12, colour = "black"))
b
#export as PDF 5X5 inches #1_N2_telomere_length_by_strand_5X5_box_violin.pdf
#export as PDF 5X3 inches #1_N2_telomere_length_by_strand_5X5_box_violin_magnified.pdf

######### plot telomere length by strand, stage and chromosome###############
e <- ggplot(df_strand_sum, aes(x = stage , y = length , fill = strand))+#by strand
  geom_violin(position = position_dodge(width = 1))+
  geom_boxplot(width = 0.2 , 
               notch = TRUE ,
               outlier.shape = NA, 
               position = position_dodge(width = 1)
  )+
  #geom_jitter(width = 0.1 , alpha = 0.1)+
  facet_grid( cols = vars(chromosome))+
  #facet_grid(chromosome~ position  )+#arrange plot by chromosome and chromosome end
  scale_y_continuous(breaks=seq(0,max_label,by = break_by),#y-axis breaks
                     labels= seq(0,max_label/1000,by = break_by/1000),#rename y-axis break labels
                     expand = c(0, 0))+#set plot to 0,0 poistion
  coord_cartesian(ylim = c(0 , max_zoom))+#show visual zoom 
  #coord_cartesian(ylim = c(2000,2500))+#0 , max_zoom))+#show visual zoom 
  #tick mark labels
  scale_x_discrete(labels=c("N2EM" = "EMB","N2L1" = "L1","N2YA" = "YA"))+
  #scale_x_discrete(labels = c("F" = "Forward","R" = "Reverse"))+
  #X-axis label
  #xlab(expression(paste(italic("C. elegans")," N2 chromosomes")))+
  ylab("Telomere sequences length (kb)")+
  #scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  #scale_fill_discrete(name = "Strand", labels = c("Forward", "Reverse"))+#legend correction
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x label
        panel.grid.minor = element_blank(),
        #panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        legend.position="top",
        legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))
e
#export as PDF 10X3.5 inches #4_N2_telomere_length_by_strand_end_chr_stage_box_violin_10X3_5_1row.pdf

######### plot telomere length by strand, position, stage and chromosome###############
d <- ggplot(df_strand_sum, aes(x = stage , y = length , fill = strand))+#by strand
  geom_violin(position = position_dodge(width = 1))+
  geom_boxplot(width = 0.2 , 
               notch = TRUE ,
               outlier.shape = NA, 
               position = position_dodge(width = 1)
  )+
  #geom_jitter(width = 0.1 , alpha = 0.1)+
  #facet_grid( cols = vars(chromosome))+
  facet_grid(chromosome~ position  )+#arrange plot by chromosome and chromosome end
  scale_y_continuous(breaks=seq(0,max_label,by = break_by),#y-axis breaks
                     labels= seq(0,max_label/1000,by = break_by/1000),#rename y-axis break labels
                     expand = c(0, 0))+#set plot to 0,0 poistion
  coord_cartesian(ylim = c(0 , max_zoom))+#show visual zoom 
  #coord_cartesian(ylim = c(2000,2500))+#0 , max_zoom))+#show visual zoom 
  #tick mark labels
  scale_x_discrete(labels=c("N2EM" = "EMB","N2L1" = "L1","N2YA" = "YA"))+
  #scale_x_discrete(labels = c("F" = "Forward","R" = "Reverse"))+
  #X-axis label
  #xlab(expression(paste(italic("C. elegans")," N2 chromosomes")))+
  ylab("Telomere sequences length (kb)")+
  #scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  scale_fill_discrete(name = "Strand", labels = c("Forward", "Reverse"))+#legend correction
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x label
        panel.grid.minor = element_blank(),
        #panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        legend.position="top",
        legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))
d
#export as PDF 10X12 inches #5_N2_telomere_length_by_strand_end_chr_position_stage_box_violin_10X12.pdf

######### plot telomere length by strand, position and chromosome###############
f <- ggplot(df_strand_sum, aes(x = chromosome , y = length , fill = strand))+#by strand
  geom_violin(position = position_dodge(width = 1))+
  geom_boxplot(width = 0.2 , 
               notch = TRUE ,
               outlier.shape = NA, 
               position = position_dodge(width = 1)
  )+
  #add mean dot to each plot
  #stat_summary(fun.y = mean, geom = "point",shape = 20,size = 2,position = position_dodge(width = 1),color = "white" ) +
  #geom_jitter(width = 0.1 , alpha = 0.1)+
  #facet_grid( cols = vars(chromosome))+
  #facet_grid(chromosome~ position  )+#arrange plot by chromosome and chromosome end
  scale_y_continuous(breaks=seq(0,max_label,by = break_by),#y-axis breaks
                     labels= seq(0,max_label/1000,by = break_by/1000),#rename y-axis break labels
                     expand = c(0, 0))+#set plot to 0,0 poistion
  coord_cartesian(ylim = c(0 , max_zoom))+#show visual zoom 
  #coord_cartesian(ylim = c(2000,2500))+#0 , max_zoom))+#show visual zoom 
  #tick mark labels
  #scale_x_discrete(labels=c("N2EM" = "EMB","N2L1" = "L1","N2YA" = "YA"))+
  #scale_x_discrete(labels = c("F" = "Forward","R" = "Reverse"))+
  #X-axis label
  #xlab(expression(paste(italic("C. elegans")," N2 chromosomes")))+
  ylab("Telomere sequences length (kb)")+
  #scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  scale_fill_discrete(name = "Strand", labels = c("Forward", "Reverse"))+#legend correction
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x label
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        legend.position="top",
        legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))
f
#export as PDF 7X4 inches #3_N2_telomere_length_by_strand_chr_box_violin_7X4_1row.pdf

#########stacked plots for reads count#############

# Stacked histogram by chromosome
g <- ggplot(df_strand_sum, aes(chromosome, fill = stage )) + 
  geom_histogram(stat="count")+#stacked histogram for counted reads by stage
  scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  scale_y_continuous(expand = c(0, 0))+
  ylab("Telomere reads count")+
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x-axis title
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        legend.position="top",
        legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))
g
#export PDF 3X5 inches
# Stacked histogram by position
aa <- ggplot(df_strand_sum, aes(chromosome, fill = stage )) + 
  geom_histogram(stat="count")+#stacked histogram for counted reads by stage
  scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  scale_y_continuous(expand = c(0, 0))+
  ylab("Telomere reads count")+
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x-axis title
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        legend.position="top",
        legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))+
  facet_grid( cols = vars(position))
aa
#export PDF 5X5 inches

# Stacked histogram by strand
ab <- ggplot(df_strand_sum, aes(chromosome, fill = stage )) + 
  geom_histogram(stat="count")+#stacked histogram for counted reads by stage
  scale_fill_discrete(name = "Stage", labels = c("EMB", "L1", "YA"))+#legend correction
  scale_y_continuous(expand = c(0, 0))+
  ylab("Telomere reads count")+
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x-axis title
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        legend.position="top",
        legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))+
  facet_grid( cols = vars(strand))
ab
#export PDF 5X5 inches

# Stacked histogram by strand
ba <- ggplot(df_strand_sum, aes(x = stage, y = length, fill = strand )) + 
  geom_bar(position="fill", stat="identity")+#stacked percent
  scale_fill_discrete(name = "Strand", labels = c("G-strand", "C-strand"))+#legend correction
  scale_x_discrete(labels=c("N2EM" = "EMB","N2L1" = "L1","N2YA" = "YA"))+#x-axis label correction
  scale_y_continuous(expand = c(0, 0),
                     breaks=seq(0,1,by = 0.25),#y-axis breaks
                     labels= seq(0,100,by = 25))+#rename y-axis break labels)
  ylab("Strand proportion %")+
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x-axis title
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        #legend.position="top",
        #legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))
ba
#PDF 3X5 inches

# Stacked histogram by strand and position
bb <- ggplot(df_strand_sum, aes(x = stage, y = length, fill = strand )) + 
  geom_bar(position="fill", stat="identity")+#stacked percent
  scale_fill_discrete(name = "Strand", labels = c("G-strand", "C-strand"))+#legend correction
  scale_x_discrete(labels=c("N2EM" = "EMB","N2L1" = "L1","N2YA" = "YA"))+#x-axis label correction
  scale_y_continuous(expand = c(0, 0),
                     breaks=seq(0,1,by = 0.25),#y-axis breaks
                     labels= seq(0,100,by = 25))+#rename y-axis break labels)
  ylab("Strand proportion")+
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x-axis title
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        #legend.position="top",
        #legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))+
  facet_grid( cols = vars(position))
bb
#PDF  4X5 inches


# Stacked histogram by chromosome
bc <- ggplot(df_strand_sum, aes(x = stage, y = length, fill = strand )) + 
  geom_bar(position="fill", stat="identity")+#stacked percent
  scale_fill_discrete(name = "Strand", labels = c("G-strand", "C-strand"))+#legend correction
  scale_x_discrete(labels=c("N2EM" = "EMB","N2L1" = "L1","N2YA" = "YA"))+#x-axis label correction
  scale_y_continuous(expand = c(0, 0),
                     breaks=seq(0,1,by = 0.25),#y-axis breaks
                     labels= seq(0,100,by = 25))+#rename y-axis break labels)
  ylab("Strand proportion")+
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),#remove x-axis title
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #legend.position = "none",#remove legend from plot
        legend.position="top",
        legend.direction = "horizontal",
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size=12, colour = "black"))+
  facet_grid( cols = vars(chromosome))
bc
#PDF  4X5 inches