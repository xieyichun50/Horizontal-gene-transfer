library(tidyr)
library(dplyr)
library(stringr)
HGT.matrix<-read.delim("Alienness_reports_v2/A)HGT_in_millipedes.txt", header = TRUE)
names(HGT.matrix)<-c("shared_species","HGTgeneid")
HGT.matrix<-separate(HGT.matrix, HGTgeneid, sep = " ", into = c("ID"), remove = FALSE)
millipedeid<-c("Ato", "Gma", "Hho", "Nno", "Tco")
for (i in 1:length(millipedeid)) {
  speciesmatrix<-read.delim(paste0("Alienness_reports_v2/", millipedeid[i], "_stat_queries_1_likely_hgt.xls"),
                            skip = 1, header = TRUE)
  speciesmatrix<-speciesmatrix[,3:4]
  names(speciesmatrix)<-c(millipedeid[i], "HGTgeneid")
  HGT.matrix<-merge(HGT.matrix, speciesmatrix, by = "HGTgeneid", all.x = TRUE)
}
write.table(HGT.matrix, "Alienness_reports_v2/A) HGT in millipedes_species.txt", 
            row.names = F, quote = F, sep = "\t")


HGT.matrix<-read.delim("Alienness_reports_v2/HGT_in_centipedes.txt", header = TRUE)
names(HGT.matrix)<-c("shared_species","HGTgeneid")
HGT.matrix<-separate(HGT.matrix, HGTgeneid, sep = " ", into = c("ID"), remove = FALSE)
centipedeid<-c("Lni", "Rim", "Sma", "Ttu")
for (i in 1:length(centipedeid)) {
  speciesmatrix<-read.delim(paste0("Alienness_reports_v2/", centipedeid[i], "_stat_queries_1_likely_hgt.xls"),
                            skip = 1, header = TRUE)
  speciesmatrix<-speciesmatrix[,3:4]
  names(speciesmatrix)<-c(centipedeid[i], "HGTgeneid")
  HGT.matrix<-merge(HGT.matrix, speciesmatrix, by = "HGTgeneid", all.x = TRUE)
}
write.table(HGT.matrix, "Alienness_reports_v2/B) HGT in centipedes_species.txt",
            row.names = F, quote = F, sep = "\t")
