library(ggplot2)
library(gggenes)
library(ggforce)

HGT<-example_genes

for (i in 1:nrow(HGT)) {
  if (HGT$orientation[i] > 0) {
    HGT$start1[i]<-HGT$start[i]
    HGT$end1[i]<-HGT$end[i]
  } else {
    HGT$start1[i]<-HGT$end[i]
    HGT$end1[i]<-HGT$start[i]
  }
}

gene.count<-as.data.frame(xtabs(~ gene, HGT))
gene.count$gene<-as.character(gene.count$gene)

for (i in 1:nrow(gene.count)) {
  if (gene.count$Freq[i] < 5) {
    gene.count$genetype[i]<-"Unique gene"
  } else {
    gene.count$genetype[i]<-gene.count$gene[i]
  }
}

HGT<-merge(HGT, gene.count, all.x = TRUE, by = "gene")

dummies <- make_alignment_dummies(
  HGT,
  aes(xmin = start, xmax = end, y = molecule, id = genetype),
  on = "genE"
)
row.names(dummies)<-1:nrow(dummies)

ggplot(HGT, 
       aes(xmin = start1, xmax = end1, y = molecule, 
           fill = genetype))+
  geom_gene_arrow()+
  geom_gene_label(aes(xmin = start1, xmax = end1, y = molecule, label = gene),
                  HGT,
                  align = "centre")+
  facet_wrap(~ molecule, scales = "free", ncol = 1)+
  geom_blank(data = dummies, mapping = aes(xmin = start, xmax = end, y = molecule))+
  labs(y = "", fill = "Genes")+
  scale_fill_brewer(palette = "Set3")+
  theme_genes()+
  theme(panel.grid.major.x = element_line(linetype = "blank"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 0, hjust = 0.5, face = "plain"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA),
        legend.position = "right")
