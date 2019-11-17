setwd("F:/nanopore_data/JUfosmid/seperation/fos_by_coordinate")
library(data.table)
for (i in 1:27){
  assign(paste0("sp",i),read.csv(paste0("cni_cutoff2/barcode",formatC(i, width = 2, flag = '0'),".csv"),sep = "\t", header = TRUE))
  assign(paste0("sp",i),na.omit(get(paste0("sp",i)),1,max))
}
for (i in 1:10){
  assign(paste0("pp",i),read.csv(paste0("cni_cutoff2/barcode",i+27,".csv"),sep = "\t", header = TRUE))
  assign(paste0("pp",i),na.omit(get(paste0("pp",i)),1,max))
}
for (i in 1:8){
  assign(paste0("rp",i),read.csv(paste0("cni_cutoff2/barcode",i+37,".csv"),sep = "\t", header = TRUE))
  assign(paste0("rp",i),na.omit(get(paste0("rp",i)),1,max))
}
for (i in 1:12){
  assign(paste0("cp",i),read.csv(paste0("cni_cutoff2/barcode",i+45,".csv"),sep = "\t", header = TRUE))
  assign(paste0("cp",i),na.omit(get(paste0("cp",i)),1,max))
}

counter = 0
hit_cutoff =100 #in range 200
#make an empty df for fosmid hit in 4 dimension library
hit_df <-data.frame(L_contig = character(),
                    L_start = integer(),
                    L_strand =character(),
                    R_contig = character(),
                    R_end = integer(),
                    R_strand = character(),
                    superpool = character(),
                    sp_confidence = integer(),
                    platepool = character(),
                    pp_confidence =integer(),
                    rowpool = character(),
                    rp_confidence = integer(),
                    columnpool = character(),
                    cp_confidence =integer(),
                    stringsAsFactors=FALSE)
for (j in 1:27){
  print(paste("Working on SuperPool", j))
  for (index in get(paste0("sp",j))$fosmid){
    sp_df<-data.frame(sp_library = character(), confidence = integer(),stringsAsFactors=FALSE)
    sp_df[nrow(sp_df)+1,] <- list(paste0("sp",j),sum(get(paste0("sp",j))[index,]$L_confidence,get(paste0("sp",j))[index,]$R_confidence))
    pp_df<-data.frame(pp_library = character(), pp_confidence = integer(),stringsAsFactors=FALSE)
    cp_df<-data.frame(cp_library = character(), cp_confidence = integer(),stringsAsFactors=FALSE)
    rp_df<-data.frame(rp_library = character(), rp_confidence = integer(),stringsAsFactors=FALSE)
    index_df <-get(paste0("sp",j))[index,]
    sp_confidence <-sum(index_df$L_confidence,index_df$R_confidence)
    hit_L_contig <- index_df$L_contig
    hit_L_location <- index_df$L_start
    hit_L_strand <- index_df$L_strand
    hit_R_contig <- index_df$R_contig
    hit_R_location <- index_df$R_end
    hit_R_strand <-index_df$R_strand
    #scanning platepool libraries
    for (p in 1:10){
      working_df <- get(paste0("pp",p))
      hit_test_df <-subset(working_df, 
                           L_contig == as.character(hit_L_contig) & 
                           L_start > hit_L_location -100 &
                           L_start < hit_L_location +100 & 
                           L_strand == as.character(hit_L_strand) &
                           R_contig == as.character(hit_R_contig) &
                           R_end > hit_R_location -100 &
                           R_end < hit_R_location +100 &
                           R_strand == as.character(hit_R_strand))
      if (nrow(hit_test_df) == 0){next}else{pp_df <-rbind(pp_df,list(pp_library=paste0("pp",p),
                                                                     pp_confidence=as.integer(sum(hit_test_df$L_confidence,
                                                                                                  hit_test_df$R_confidence))),
                                                          stringsAsFactors=FALSE)}}
    for (r in 1:8){
      working_df <- get(paste0("rp",r))
      hit_test_df <-subset(working_df, 
                           L_contig == as.character(hit_L_contig) & 
                             L_start > hit_L_location -100 &
                             L_start < hit_L_location +100 & 
                             L_strand == as.character(hit_L_strand) &
                             R_contig == as.character(hit_R_contig) &
                             R_end > hit_R_location -100 &
                             R_end < hit_R_location +100 &
                             R_strand == as.character(hit_R_strand))
      if (nrow(hit_test_df) == 0){next}else{rp_df <-rbind(rp_df,list(rp_library=paste0("rp",r),
                                                                     rp_confidence=as.integer(sum(hit_test_df$L_confidence,
                                                                                                  hit_test_df$R_confidence))),
                                                                     stringsAsFactors=FALSE)}}
    for (c in 1:12){
      working_df <- get(paste0("cp",c))
      hit_test_df <-subset(working_df, 
                           L_contig == as.character(hit_L_contig) & 
                             L_start > hit_L_location -100 &
                             L_start < hit_L_location +100 & 
                             L_strand == as.character(hit_L_strand) &
                             R_contig == as.character(hit_R_contig) &
                             R_end > hit_R_location -100 &
                             R_end < hit_R_location +100 &
                             R_strand == as.character(hit_R_strand))
      if (nrow(hit_test_df) == 0){next}else{cp_df <-rbind(cp_df,list(cp_library=paste0("cp",c),
                                                                     cp_confidence=as.integer(sum(hit_test_df$L_confidence,
                                                                                                  hit_test_df$R_confidence))),
                                                                     stringsAsFactors=FALSE)}}
    if (nrow(pp_df) ==0 |nrow(rp_df)==0 |nrow(cp_df)==0){
      next
    }else{
      for (pp_element in pp_df$pp_library) {
        for(rp_element in rp_df$rp_library){
          for (cp_element in cp_df$cp_library){
            counter = counter +1
            hit_df<-rbind(hit_df,
                          list(L_contig = as.character(index_df$L_contig),
                               L_start = index_df$L_start,
                               L_strand = as.character(index_df$L_strand),
                               R_contig = as.character(index_df$R_contig),
                               R_end = index_df$R_end,
                               R_strand = as.character(index_df$R_strand),
                               superpool = paste0("sp",j),
                               sp_confidence = sp_confidence,
                               platepool = pp_element,
                               pp_confidence = as.numeric(as.character(pp_df[pp_df$pp_library==pp_element,]$pp_confidence)),
                               rowpool = rp_element,
                               rp_confidence = as.numeric(as.character(rp_df[rp_df$rp_library==rp_element,]$rp_confidence)),
                               columnpool = cp_element,
                               cp_confidence = as.numeric(as.character(cp_df[cp_df$cp_library==cp_element,]$cp_confidence))),
                          stringsAsFactors=FALSE)
          }
        }
      }
    }
  }
}
print(paste("Find ",counter, "in all libraries."))
write.table(hit_df,file=paste0('cni_hit_cutoff2.csv'),
            col.names=T,row.names = F,sep='\t')



