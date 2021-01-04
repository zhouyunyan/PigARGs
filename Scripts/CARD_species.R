######Distribution of ARGs in different pig farms############
data <- read.table(file = "species24_gene194_abun.xls",header = T,sep = "\t",check.names=F)

#transform to 0,1
for (i in 7:ncol(data)) {
  col <- data[,i]
  data[which(col >0),i] <- 1
  data[which(col==0),i] <- 0
}

#another method:use vegan package 
#library(vegan)
#table <- decostand(data[,7:200], method="pa")

library(pheatmap)
library(RColorBrewer)
data$Farm <- factor(data$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))

#add annotation to 194 genes
species_gene <- read.delim("clipboard",header = T,row.names = 1, check.names = F)#CARD_species.xlsx:gene_species_194
species <- as.data.frame(species_gene$species)
rownames(species) <- rownames(species_gene)
names(species) <- "Host of genes"

data.sort <- data[order(data$Farm),]
#add annotation to 425 samples
Farm <- as.data.frame(data.sort[,6])
rownames(Farm) <- data.sort$ID
names(Farm) <- "Farm"

data.pheatmap <- as.data.frame(t(data.sort[,c(7:200)]),check.names=F)
colnames(data.pheatmap) <- data.sort$ID
data.pheatmap <- data.pheatmap[match(rownames(species),rownames(data.pheatmap)),]
p1 <- pheatmap(data.pheatmap,
         cluster_rows =F,cluster_cols =F,legend = F,
         show_rownames=F,show_colnames=F,
         annotation_row = species,annotation_col = Farm)


tiff(filename = "gene194.Farm.pheatmap1.tif",width = 4800,height = 3600,res=600,compression="lzw")
p
dev.off()

#gene cluster
p2 <- pheatmap(data.pheatmap,
               cluster_rows =T,cluster_cols =F,legend = F,
               show_rownames=F,show_colnames=F,
               annotation_row = species,annotation_col = Farm)


tiff(filename = "gene194.Farm.pheatmap.cluster1.tif",width = 4800,height = 4800,res=600,compression="lzw")
p2
dev.off()

##The genes of different species were ploted separately
data.pheatmap$GeneID <- rownames(data.pheatmap)
species_gene <- read.delim("clipboard",header = T, check.names = F)#CARD_species.xlsx:gene_species_194
data.add.Species <- merge(species_gene,data.pheatmap,by="GeneID")
annotation_row <- data.add.Species[which(data.add.Species$species=="Escherichia_coli"),c(2,3)]

##Escherichia_coli
data.Escherichia_coli <- data.add.Species[which(data.add.Species$species=="Escherichia_coli"),6:430]

p.E.coli <- pheatmap(data.Escherichia_coli,
         cluster_rows =T,cluster_cols =F,legend = F,
         show_rownames=F,show_colnames=F,
         annotation_col = Farm)
tiff(filename = "E.coli.pheatmap.tif",width = 4800,height = 3600,res=600,compression="lzw")
p.E.coli
dev.off()

##Acinetobacter_baumannii
data.Acinetobacter_baumannii <- data.add.Species[which(data.add.Species$species=="Acinetobacter_baumannii"),6:430]
p.A.baumannii <- pheatmap(data.Acinetobacter_baumannii,
                     cluster_rows =T,cluster_cols =F,legend = F,
                     show_rownames=F,show_colnames=F,
                     annotation_col = Farm)
tiff(filename = "A.baumannii.pheatmap.tif",width = 4800,height = 1000,res=600,compression="lzw")
p.A.baumannii
dev.off()

##Staphylococcus_aureus
data.Staphylococcus_aureus <- data.add.Species[which(data.add.Species$species=="Staphylococcus_aureus"),6:430]
p.Staphylococcus_aureus <- pheatmap(data.Staphylococcus_aureus,
                          cluster_rows =T,cluster_cols =F,legend = F,
                          show_rownames=F,show_colnames=F,
                          annotation_col = Farm)
tiff(filename = "Staphylococcus_aureus.pheatmap.tif",width = 4800,height = 1200,res=600,compression="lzw")
p.Staphylococcus_aureus
dev.off()

#Klebsiella_pneumoniae
data.Klebsiella_pneumoniae <- data.add.Species[which(data.add.Species$species=="Klebsiella_pneumoniae"),6:430]
p.Klebsiella_pneumoniae <- pheatmap(data.Klebsiella_pneumoniae,
                                    cluster_rows =T,cluster_cols =F,legend = F,
                                    show_rownames=F,show_colnames=F,
                                    annotation_col = Farm)
tiff(filename = "Klebsiella_pneumoniae.pheatmap.tif",width = 4800,height = 600,res=600,compression="lzw")
p.Klebsiella_pneumoniae
dev.off()

#Enterococcus_faecium
data.Enterococcus_faecium <- data.add.Species[which(data.add.Species$species=="Enterococcus_faecium"),6:430]
p.Enterococcus_faecium <- pheatmap(data.Enterococcus_faecium,
                                   cluster_rows =T,cluster_cols =F,legend = F,
                                   show_rownames=F,show_colnames=F,
                                   annotation_col = Farm)
tiff(filename = "Enterococcus_faecium.pheatmap.tif",width = 4800,height = 500,res=600,compression="lzw")
p.Enterococcus_faecium
dev.off()

#Enterococcus_faecalis
data.Enterococcus_faecalis <- data.add.Species[which(data.add.Species$species=="Enterococcus_faecalis"),6:430]
p.Enterococcus_faecalis <- pheatmap(data.Enterococcus_faecalis,
                                    cluster_rows =T,cluster_cols =F,legend = F,
                                    show_rownames=F,show_colnames=F,
                                    annotation_col = Farm)
tiff(filename = "Enterococcus_faecalis.pheatmap.tif",width = 4800,height = 400,res=600,compression="lzw")
p.Enterococcus_faecalis
dev.off()

#Probiotics
data.Probiotics <- data.add.Species[which(data.add.Species$species %in% c("Akkermansia_muciniphila","Bacillus_subtilis","Bacteroides_fragilis","Faecalibacterium_prausnitzii","Lactococcus_lactis")),6:430]
p.Probiotics <- pheatmap(data.Probiotics,
                         cluster_rows =T,cluster_cols =F,legend = F,
                         show_rownames=F,show_colnames=F,
                         annotation_col = Farm)
tiff(filename = "Probiotics.pheatmap.tif",width = 4800,height = 500,res=600,compression="lzw")
p.Probiotics
dev.off()

#output gene name order by showing in the heatmap
E.coli.label <- rownames(data.Escherichia_coli)[ p.E.coli$tree_row$order]
A.baumannii.label <- rownames(data.Acinetobacter_baumannii)[ p.A.baumannii$tree_row$order]
Staphylococcus_aureus.label <- rownames(data.Staphylococcus_aureus)[ p.Staphylococcus_aureus$tree_row$order]
Klebsiella_pneumoniae.label <- rownames(data.Klebsiella_pneumoniae)[ p.Klebsiella_pneumoniae$tree_row$order]
Enterococcus_faecium.label <- rownames(data.Enterococcus_faecium)[ p.Enterococcus_faecium$tree_row$order]
Enterococcus_faecalis.label <- rownames(data.Enterococcus_faecalis)[ p.Enterococcus_faecalis$tree_row$order]
Probiotics.label <- rownames(data.Probiotics)[ p.Probiotics$tree_row$order]

labels_sort <- c(E.coli.label,A.baumannii.label,Staphylococcus_aureus.label,Klebsiella_pneumoniae.label,
                 Enterococcus_faecium.label,Enterococcus_faecalis.label,Probiotics.label)
labels_sort.table <- data.add.Species[match(labels_sort,rownames(data.add.Species)),1:5]
write.csv(labels_sort.table,file = "pheatmap.cluster.label.csv")
