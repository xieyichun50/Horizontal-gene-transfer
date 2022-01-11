library(tidyr)
library(dplyr)
library(stringr)

gtf.file<-as.data.frame(matrix(NA, nrow = 1, ncol = 5))
names(gtf.file)<-c("chr", "start", "end", "strand", "gene")

setwd("D:/centipede/HGT/seqids/centipede")

HGT.matrix<-read.delim("D:/centipede/HGT/Alienness_reports_v2/HGT_in_centipedes.txt", header = TRUE)
names(HGT.matrix)<-c("shared_species","HGTgeneid")
HGT.matrix$ID<-HGT.matrix$HGTgeneid
HGT.matrix$ID<-gsub(" / none", "", HGT.matrix$ID)
HGT.matrix$ID<-gsub(" / ", "_", HGT.matrix$ID)
HGT.matrix$ID<-gsub(" ", "", HGT.matrix$ID)
#HGT.matrix<-separate(HGT.matrix, HGTgeneid, sep = " / ", into = c("ID","ID2"), remove = FALSE)
write.table(HGT.matrix, "HGT.matrix.centipede.txt", sep = "\t", row.names = F, quote = F)

species<-read.delim("centipedes.txt", header = TRUE)

setwd("D:/centipede/HGT/seqids/millipede")
HGT.matrix<-read.delim("D:/centipede/HGT/Alienness_reports_v2/HGT_in_millipedes.txt", header = TRUE)
names(HGT.matrix)<-c("shared_species","HGTgeneid")
HGT.matrix$ID<-HGT.matrix$HGTgeneid
HGT.matrix$ID<-gsub(" / none", "", HGT.matrix$ID)
HGT.matrix$ID<-gsub(" / ", "_", HGT.matrix$ID)
HGT.matrix$ID<-gsub(" ", "", HGT.matrix$ID)
#HGT.matrix<-separate(HGT.matrix, HGTgeneid, sep = " / ", into = c("ID","ID2"), remove = FALSE)
write.table(HGT.matrix, "HGT.matrix.millipedes.txt", sep = "\t", row.names = F, quote = F)

species<-read.delim("millipedes.txt", header = TRUE)

for (i in 1:nrow(HGT.matrix)) {
  output.dir=paste0(HGT.matrix$shared_species[i], "/", HGT.matrix$ID[i])
  if (!dir.exists(output.dir)){
    dir.create(file.path(output.dir), recursive =T)
  } else {
    print("Dir already exists!")
  }
  
  for (j in 1:nrow(species)) {
    HGTgene<-read.delim(paste0("D:/centipede/HGT/Alienness_reports_v2/", 
                               species$id[j], "_stat_queries_1_likely_hgt.xls"),
                        skip = 1, header = TRUE)
    HGTgene<-HGTgene[HGTgene$best_donor_acc...best_toi_acc == HGT.matrix$HGTgeneid[i],]
    if (nrow(HGTgene) > 0) {
      HGTgene<-HGTgene[,3:4]
      
      gtf<-read.delim(paste0("D:/centipede/gtf/", species$species[j], ".longest-gene.gtf.id.gff3"),
                      header = FALSE)
      gtf$V9<-gsub("ID=","", gtf$V9)
      gtf<-separate(gtf, V9, into = c("ID"), sep = ";", remove = TRUE)
      gtf<-gtf[gtf$V3 == "mRNA", c(1,4,5,7,9)]
      names(gtf)<-c("chr", "start", "end", "strand", "gene")
      
      gtf.species<-gtf.file
      
      for (k in 1:nrow(HGTgene)) {
        gtf.gene<-gtf[gtf$gene == HGTgene$query[k],]
        gtf.sub<-gtf[gtf$chr == gtf.gene$chr,]
        
        if (nrow(gtf.sub)>1) {
          gtf.sub<-gtf.sub[order(gtf.sub$start),]
          row.names(gtf.sub)<-1:nrow(gtf.sub)
          line<-as.numeric(row.names(gtf.sub[which(gtf.sub$gene == HGTgene$query[k]),]))
          if (line > 5){
            gtf.sub1<-gtf.sub[(line-5):(line+5), ]
            
          } else {
            gtf.sub1<-gtf.sub[1:(line+5), ]
          }
          
        } else {
          gtf.sub1<-gtf.sub
        }
        
        gtf.species<-rbind(gtf.species, gtf.sub1)
      }
      gtf.species<-gtf.species[is.na(gtf.species$chr) == F, ]
      
      if (nrow(gtf.species)>0) {
        gtf.species<-unique(gtf.species)
        gtf.species$species<-species$species[j]
        gtf.species$genetype[gtf.species$gene %in% HGTgene$query]="HGT"
        write.table(gtf.species, 
                    paste0(output.dir,"/", species$species[j]),
                    sep = "\t", row.names = F, quote = F)
      } else {
        
      }
      
    } else {
      
    }
    
  }
    
}

