library(tidyr)
library(dplyr)
library(stringr)

gtf.file<-as.data.frame(matrix(NA, nrow = 1, ncol = 5))
names(gtf.file)<-c("chr", "start", "end", "strand", "gene")

geneno=5

####centipede
setwd("/jelly_data/yichun/myriapod/HGT/seqids/centipede")
HGT.matrix<-read.delim("HGT.matrix.centipedes.txt", header = TRUE)
species<-read.delim("centipedes.txt", header = TRUE)

for (i in 1:nrow(HGT.matrix)) {
  output.dir=paste0(HGT.matrix$shared_species[i], "/", HGT.matrix$ID[i])
  if (!dir.exists(output.dir)){
    dir.create(file.path(output.dir), recursive =T)
  } else {
    print("Dir already exists!")
  }
  
  for (j in 1:nrow(species)) {
    HGTgene<-read.delim(paste0("/jelly_data/yichun/myriapod/HGT/alienness/", 
                               species$id[j], "_stat_queries_1_likely_hgt.xls"),
                        skip = 1, header = TRUE)
    HGTgene<-HGTgene[HGTgene$best_donor_acc...best_toi_acc == HGT.matrix$HGTgeneid[i],]
    if (nrow(HGTgene) > 0) {
      HGTgene<-HGTgene[,3:4]
      
      gtf<-read.delim(paste0("/data/huilab_genome/centipede_genome/gene_modele_final/", species$species[j], ".longest-gene.gtf.id.gff3"),
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
          if (line > geneno){
            gtf.sub1<-gtf.sub[(line-geneno):(line+geneno), ]
            
          } else {
            gtf.sub1<-gtf.sub[1:(line+geneno), ]
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

####millipede
setwd("/jelly_data/yichun/myriapod/HGT/seqids/millipede")
HGT.matrix<-read.delim("HGT.matrix.millipedes.txt", header = TRUE)
species<-read.delim("millipedes.txt", header = TRUE)

for (i in 1:nrow(HGT.matrix)) {
  output.dir=paste0(HGT.matrix$shared_species[i], "/", HGT.matrix$ID[i])
  if (!dir.exists(output.dir)){
    dir.create(file.path(output.dir), recursive =T)
  } else {
    print("Dir already exists!")
  }
  
  for (j in 1:nrow(species)) {
    HGTgene<-read.delim(paste0("/jelly_data/yichun/myriapod/HGT/alienness/", 
                               species$id[j], "_stat_queries_1_likely_hgt.xls"),
                        skip = 1, header = TRUE)
    HGTgene<-HGTgene[HGTgene$best_donor_acc...best_toi_acc == HGT.matrix$HGTgeneid[i],]
    if (nrow(HGTgene) > 0) {
      HGTgene<-HGTgene[,3:4]
      
      gtf<-read.delim(paste0("/data/huilab_genome/centipede_genome/gene_modele_final/", species$species[j], ".longest-gene.gtf.id.gff3"),
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
          if (line > geneno){
            gtf.sub1<-gtf.sub[(line-geneno):(line+geneno), ]
            
          } else {
            gtf.sub1<-gtf.sub[1:(line+geneno), ]
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

