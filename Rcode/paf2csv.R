setwd("F:/nanopore_data/JUfosmid/seperation")
fileseq <- seq(1,24)
header = c('qname', 'qlen','qstart','qend','strand','rname','rlen','rstart','rend','match','alignlength')
for (i in fileseq){
  fileName = paste0('align/trim_barcode',formatC(i, width = 2, flag = '0'),'.paf')
  saving = paste0('paf2csv/filter_barcode',formatC(i, width = 2, flag = '0'),'.csv')
  df <- read.csv(file = fileName, header = F, sep = '\t')
  df_ex <-subset(df, select = c(seq(1,11)))
  #in some cases, left flanking or right flanking is long enough to align to other flanking sequence
  df_ex2 <- subset(df_ex, df_ex$V1 != paste0(df_ex$V6))
  colnames(df_ex2) <-header
  df_ex2$qper <- with(df_ex2, abs(qend-qstart)/qlen )
  df_ex2$rper <- with(df_ex2, abs(rend-rstart)/rlen )
  write.table(df_ex2, file = saving, sep = '\t',row.names=FALSE)
  print(paste("Complete", formatC(i , width = 2 , flag = '0')))
}
##############function test####################
library(readr)
library(data.table)
filter_barcode03 <- read_delim("paf2csv/filter_barcode03.csv","\t", escape_double = FALSE, trim_ws = TRUE)

sample <- subset(filter_barcode03, strand == '+' & (qper >0.80 | rper >0.80))
read_id = '44b73770-d1c6-40f4-979b-8e069e21c2a8'
tmp1 <- sample[sample$qname %like% read_id, ]
tmp2 <- sample[sample$rname %like% read_id, ]
tmp_merge <- rbind(tmp1,tmp2)
View(tmp_merge)

tmp1 <- sample[sample$qname %like% "146b2a04-7ce6-4d7b-8d45-8d4e2adc259d", ]
tmp2 <- sample[sample$rname %like% "146b2a04-7ce6-4d7b-8d45-8d4e2adc259d", ]
tmp3 <- sample[sample$qname %like% "ddc38594-7e2f-415a-9c04-a8addfeafef7", ]
tmp4 <- sample[sample$qname %like% "ddc38594-7e2f-415a-9c04-a8addfeafef7", ]
merge <-rbind(tmp1,tmp2,tmp3,tmp4)
View(merge)
#convert column in alignment table to uniq list
uniq_list_1 = c()
uniq_list_2 = c()
for (readName in CJ(merge$qname,  unique = T)){uniq_list_1 <- readName}
for (readName in CJ(merge$rname,  unique = T)){uniq_list_2 <- readName}
uniq_list<-c(uniq_list_1,uniq_list_2)
View(uniq_list)

left_list <-grep( pattern = "L",uniq_list , value = T)
right_list <-grep( pattern = "R",uniq_list, value = T)
View(left_list)
View(right_list)
tmpLst<- list(left_list,right_list)
names(tmpLst) <- c("left", "right")
View(tmpLst)
class(tmpLst)

fos_counter = 0
j = 3

# read in barcoded reads with backbone sequence, extract their read name
fileName = paste0('trim/trim_barcode',formatC(j, width = 2, flag = '0'),'.list')
readInFile = file(fileName, "r")
readList <-readInFile
read_num = 0
while (TRUE) {
  line = readLines(readList, n=1)
  if(length(line)==0){break}
  read_name = paste0(line)
  read_num = read_num + 1
}






readlist_file <- read.csv(file = fileName, header = F, sep = '\t')
readlist = c()
for (element in CJ(readlist_file$V1,  unique = T)){readlist <- element}

for (j in seq(1, length(readlist_file$V1))){
  readname = readlist_file[j,]
  print(readname)
}





#make a empty df
fosmid_df <-data.frame(counter = integer(),
                       readname = factor())
fosmid_df



uniq_list

uniq_list2
uniq_list3 <-list(uniq_list)

fos_list

names(uniq_list3) <- paste0("fos",fosmid_counter+1)
fos_list <- c(fos_list,uniq_list3)
View(fos_list)



########test substring#######

View(uniq_list_tmp)
sub_test<- substr(uniq_list_tmp,start = 1 , stop = 36)
View(sub_test)
sub_test
uniq_sub_test<-unique(sub_test)
View(uniq_sub_test)


View(fos_list)

i=1
for (list in fos_list){
  print(paste('length:',length(list)))
  i=i+1
  if (i >6){break}
}

fin_fos_list = list()
cutoff = 3
R_counter = 1
for (list in fos_list){
  if(length(list) < cutoff){
    next
  }else{
    shrink_list <- list(list)
    names(shrink_list) <- paste0("R_fos",R_counter)
    fin_fos_list <-c(fin_fos_list , shrink_list)
    R_counter = R_counter + 1
  }
}

View(fin_fos_list)

require(seqinr)

test_fa <-read.fasta("trim/trim_barcode03.fa", as.string = T)
View(test_fa)
class(test_fa)

read_name = 'c5c089c6-8f65-41ba-a0d1-145ae651ac55'

tmp1 <- readInFilter[readInFilter$qname %like% sample_name, ]
tmp2 <- readInFilter[readInFilter$rname %like% sample_name, ]

grep_sample <- rbind(tmp1,tmp2)
View(grep_sample)
