#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gggenes))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(RColorBrewer))

##Create parameters
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="input directory' [default %default]",
              dest="input.dir"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="output filename dir (end with /) [default %default]",
              dest="output.dir"),
  make_option(c("-s","--share"), type="logical", default=FALSE,
              help="show scaffolds with orthogroups only [default %default]",
              dest="share"),
  make_option(c("-p","--paralog"), type="logical", default=FALSE,
              help="show paralogs [default %default]",
              dest="paralog")
)
options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

#opt$input.dir="D:/centipede/HGT/seqids/millipede/Ato_Nno/WP_029904946.1"
#opt$share=TRUE
#opt$input.dir="D:/centipede/HGT/seqids/centipede/Lni_Rim_Sma_Ttu/WP_189008864.1"
#opt$output.dir="D:/centipede/HGT/seqids/figure/"

##filename prefix
outfilename<-opt$input.dir[1]
outfilename<-as.data.frame(strsplit(outfilename, "/"))
opt$group<-outfilename[nrow(outfilename)-2,1]
opt$species<-outfilename[nrow(outfilename)-1,1]
opt$HGT<-outfilename[nrow(outfilename),1]
opt$outfilename<-paste0(opt$group,"_",opt$species,"_",opt$HGT)
rm(outfilename)

##species name table
if (opt$group == "centipede") {
  #species<-read.delim("D:/centipede/HGT/seqids/centipede/centipedes.txt", header = TRUE)
  species<-read.delim("/jelly_data/yichun/myriapod/HGT/seqids/centipede/centipedes.txt", header = TRUE)
} else if (opt$group == "millipede") {
  #species<-read.delim("D:/centipede/HGT/seqids/millipede/millipedes.txt", header = TRUE)
  species<-read.delim("/jelly_data/yichun/myriapod/HGT/seqids/millipede/millipedes.txt", header = TRUE)
} else if (opt$group == "common"){
  #species<-read.delim("D:/centipede/HGT/seqids/common/common.txt", header = TRUE)
  species<-read.delim("/jelly_data/yichun/myriapod/HGT/seqids/common/common.txt", header =TRUE)
} else {}

specieslist<-opt$species[1]
specieslist<-as.data.frame(strsplit(specieslist, "_"))
names(specieslist)<-c("id")
specieslist<-merge(specieslist, species, by = "id", all.x = TRUE)

##location table
HGT<-as.data.frame(matrix(NA, nrow = 1, ncol = 7))
names(HGT)<-c("chr", "start", "end", "strand", "gene", "species", "genetype")

for (i in 1:nrow(specieslist)) {
  HGT.species<-NA
  HGT.species<-read.delim(paste0(opt$input.dir[1], "/", specieslist$species[i]))
  reversechr<-unique(HGT.species[HGT.species$genetype=="HGT" & HGT.species$strand == "-" , "chr"])
  HGT.species[HGT.species$chr %in% reversechr, "start"]<- 0-HGT.species[HGT.species$chr %in% reversechr, "start"]
  HGT.species[HGT.species$chr %in% reversechr, "end"]<- 0-HGT.species[HGT.species$chr %in% reversechr, "end"]
  HGT<-rbind(HGT, HGT.species)
}

HGT<-HGT[is.na(HGT$chr)==F,]
HGT<-unique(HGT)

#Orthogroup ID
OG<-read.delim(paste0(opt$input.dir[1], "/Orthogroups.txt"), header = F)
OGcol<-paste0("gene", 1:nrow(as.data.frame(strsplit(OG$V1[1], " "))))

OG<-separate(OG, "V1", sep = " ", into = c("OG", OGcol))
OG$OGname<-NA
if (opt$paralog == F) {
  for (i in 1:nrow(OG)) {
    nsp<-substr(as.list(OG[i,2:(ncol(OG)-1)]), 1, 4)
    nsp<-nsp[is.na(nsp)==F]
    if (length(unique(nsp))<=1) {
      OG$OGname[i]="Unique gene"
    } else {
      OG[is.na(OG$gene2)==T, "OGname"]="Unique gene"
    }
  }
}
OG<-OG[order(OG$OGname, decreasing = F, na.last = F),]
row.names(OG)<-1:nrow(OG)
OG$OGname[is.na(OG$OGname)]=paste0("Orthogroup", LETTERS[as.numeric(row.names(OG[is.na(OG$OGname)==T,]))])

OG<-OG[,2:ncol(OG)]

OGgene<-as.data.frame(matrix(NA, nrow = 1, ncol = 2))
names(OGgene)=c("OGname", "gene")
for (i in 1:nrow(OG)) {
  OGgene.sub<-NA
  OGgene.sub<-as.data.frame(t(OG[i,1:(ncol(OG)-2)]))
  names(OGgene.sub)[1]="gene"
  OGgene.sub$OGname<-OG$OGname[i]
  OGgene<-rbind(OGgene, OGgene.sub)
}
rm(OGgene.sub)
OGgene<-OGgene[is.na(OGgene$gene)==F,]
OGgene<-unique(OGgene)
row.names(OGgene)<-1:nrow(OGgene)

##Merge genetype information
HGT<-merge(HGT, OGgene, by = "gene", all = TRUE)
HGT[is.na(HGT$genetype)==F, "genetype"]<-paste0(HGT[is.na(HGT$genetype)==F, "genetype"],"(",opt$HGT[1],")")
HGT[is.na(HGT$genetype), "genetype"]<-HGT[is.na(HGT$genetype), "OGname"]

for (i in 1:nrow(HGT)) {
  if (HGT$strand[i] == "+") {
    HGT$start1[i]<-HGT$start[i]
    HGT$end1[i]<-HGT$end[i]
  } else {
    HGT$start1[i]<-HGT$end[i]
    HGT$end1[i]<-HGT$start[i]
  }
}

##Get scaffolds with orthologs
HGT$specieschr<-paste0(HGT$species,"\n",HGT$chr)
if (opt$share == TRUE) {
  chrsub<-unique(HGT[HGT$genetype != paste0("HGT","(",opt$HGT[1],")") & 
                       HGT$genetype != "Unique gene", "specieschr"])
  HGT<-HGT[HGT$specieschr %in% chrsub,]
} else {
  
}

##plot
if (nrow(HGT)>0) {
  dummies <- make_alignment_dummies(HGT,
                                    aes(xmin = start, xmax = end, y = specieschr, id = genetype),
                                    on = paste0("HGT(",opt$HGT[1],")"))
  row.names(dummies)<-1:nrow(dummies)
  dummies<-dummies[!duplicated(dummies$specieschr),]
  
  ##set xscale
  yscale<-ceiling(log10(max(abs(dummies[,2:3]))))
  
  if (yscale > 6) {
    yscale<-10^6
    xlab = "(MB)"
  } else {
    yscale<-10^3
    xlab = "(KB)"
  }
  
  ##setcolor
  colorset3<-brewer.pal(12, "Set3")
  if (length(unique(HGT$genetype))<=2) {
    genecolor=c("Red","white")
  } else {
    colorset<-colorset3[1:(length(unique(HGT$genetype))-2)]
    genecolor=c("Red", colorset, "white")
  }
  
  write.table(HGT,
              paste0(opt$output.dir,opt$outfilename,".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  ggplot(HGT, 
         aes(xmin = start1/yscale, xmax = end1/yscale, y = specieschr, 
             fill = genetype))+
    geom_gene_arrow()+
    geom_gene_label(aes(xmin = start1/yscale, xmax = end1/yscale, y = specieschr, label = gene),
                    HGT,
                    align = "centre")+
    facet_wrap(~ specieschr, scales = "free", ncol = 1)+
    geom_blank(data = dummies, mapping = aes(xmin = start/yscale, xmax = end/yscale, y = specieschr))+
    labs(x= xlab, y = "", fill = "Gene groups")+
    scale_fill_manual(values = genecolor)+
    theme_genes()+
    theme(panel.grid.major.x = element_line(linetype = "blank"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          #plot.title = element_text(size = 12, hjust = 0.5, face = "plain"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = "right")
  
  ggsave(paste0(opt$output.dir,opt$outfilename,"_scaled.png"),
         width = 20, height = 1*nrow(dummies), units = "in", dpi = 300, limitsize = F)
  
  ggplot(HGT, 
         aes(xmin = start1/yscale, xmax = end1/yscale, y = specieschr, 
             fill = genetype))+
    geom_gene_arrow()+
    geom_gene_label(aes(xmin = start1/yscale, xmax = end1/yscale, y = specieschr, label = gene),
                    HGT,
                    align = "centre")+
    facet_wrap(~ specieschr, scales = "free", ncol = 1)+
    #geom_blank(data = dummies, mapping = aes(xmin = start/yscale, xmax = end/yscale, y = specieschr))+
    labs(x= xlab, y = "", fill = "Gene groups")+
    scale_fill_manual(values = genecolor)+
    theme_genes()+
    theme(panel.grid.major.x = element_line(linetype = "blank"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12, colour = "black"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          #plot.title = element_text(size = 12, hjust = 0.5, face = "plain"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key = element_rect(fill = NA),
          legend.background = element_rect(fill = NA),
          legend.position = "right")
  
  ggsave(paste0(opt$output.dir,opt$outfilename,"_unscale.png"),
         width = 20, height = 1*nrow(dummies), units = "in", dpi = 300, limitsize = F)
} else {
  cat(paste0(opt$output.dir,opt$outfilename, " is empty!"))
}


