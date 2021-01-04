#import data£ºARO.total.xlsx
Freq <- read.delim("clipboard",header = T,check.names = F) ##sheet:ARO

######### Number of ARGs ###############
table(Freq$Farm)
# Dingnan  Jiangying    KD-1400    KD-3480    KD-3800      NC-F6 NC-Tibetan 
# 63         16          8          6          7        293          6 
# Shahu       Wild
# 20          6 

Freq$Farm <- factor(Freq$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))
library(ggpubr)
#pairwise comparison
compare <- compare_means(ARO_Freq~Farm, data=Freq,p.adjust.method = "fdr")
#multiple group comparison
compare_means(ARO_Freq~Farm, data=Freq,p.adjust.method = "fdr",method = "kruskal.test")
# .y.             p    p.adj p.format p.signif method        
# <chr>       <dbl>    <dbl> <chr>    <chr>    <chr>         
#   1 ARO_Freq 5.56e-24 5.56e-24 <2e-16   ****     Kruskal-Wallis

tiff(filename = "Farm.ARO.Freq.tif",width = 2500,height = 1800,res=600,compression="lzw")
ggboxplot(Freq, x="Farm", y="ARO_Freq",  fill  = "Farm",palette = "npg") + 
  guides(fill=FALSE)+
  labs(y="AMR gene numbers",x=NULL) + 
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title = element_text(size=10))
dev.off()

######### Abundance of ARGs ###############
compare_means(Abun~Farm, data=Freq,p.adjust.method = "fdr",method = "kruskal.test")
# .y.          p    p.adj p.format p.signif method        
# <chr>    <dbl>    <dbl> <chr>    <chr>    <chr>         
#   1 Abun  1.70e-18 1.70e-18 <2e-16   ****     Kruskal-Wallis

tiff(filename = "Farm.ARO.Abun.tif",width = 2500,height = 1800,res=600,compression="lzw")
ggboxplot(Freq, x="Farm", y="Abun",  fill  = "Farm",palette = "npg") + 
  guides(fill=FALSE)+
  labs(y="AMR abundance (FPKM)",x=NULL) + 
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title = element_text(size=10))
dev.off()

######### Shannon ###############
compare_means(species_Shannon~Farm, data=Freq,p.adjust.method = "fdr",method = "kruskal.test")
# .y.                    p    p.adj p.format p.signif method        
# <chr>              <dbl>    <dbl> <chr>    <chr>    <chr>         
#   1 species_Shannon 6.74e-21 6.74e-21 <2e-16   ****     Kruskal-Wallis

tiff(filename = "Farm.species.Shannon.tif",width = 2500,height = 1800,res=600,compression="lzw")
ggboxplot(Freq, x="Farm", y="species_Shannon",  fill  = "Farm",palette = "npg") + 
  guides(fill=FALSE)+
  labs(y="Shannon index (species level)",x=NULL) + 
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title = element_text(size=10))
dev.off()

######### species number###############
compare_means(species_number~Farm, data=Freq,p.adjust.method = "fdr",method = "kruskal.test")
# .y.                   p    p.adj p.format p.signif method        
# <chr>             <dbl>    <dbl> <chr>    <chr>    <chr>         
#   1 species_number 2.14e-43 2.14e-43 <2e-16   ****     Kruskal-Wallis
tiff(filename = "Farm.species.number.tif",width = 2500,height = 1800,res=600,compression="lzw")
ggboxplot(Freq, x="Farm", y="species_number",  fill  = "Farm",palette = "npg") + 
  guides(fill=FALSE)+
  labs(y="Number of species (species level)",x=NULL) + 
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title = element_text(size=10))
dev.off()

######### Simpson ###############
tiff(filename = "Farm.species.Simpson.tif",width = 2500,height = 1800,res=600,compression="lzw")
ggboxplot(Freq, x="Farm", y="species_Simpson",  fill  = "Farm",palette = "npg") + 
  guides(fill=FALSE)+
  labs(y="Simpson's index (species level)",x=NULL) + 
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.title = element_text(size=10))
dev.off()


#########Beta diversity###############
data <- read.table("adult425.faeces.ARO.xls",header = T,sep = "\t",check.names = F)
data.Abun <- data[,-c(1:6)]
rownames(data.Abun) <-data$ID
library(vegan)
data.Abun.hel <- decostand(data.Abun, "hellinger")

distance <- vegdist(data.Abun.hel,method="bray") 
pcoa = cmdscale(distance, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("PCoA1", "PCoA2", "PCoA3") 
eig = pcoa$eig
points = cbind(points,data$Farm)
names(points)[4] <- "Farm" 

points$Farm <- factor(points$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))


library(ggpubr)
p <- ggscatter(points,x="PCoA1", y="PCoA2", size = 0.8,
               color="Farm",ellipse = TRUE,  
               mean.point = TRUE, star.plot = TRUE,   
               ellipse.level = 0.9,  
               ggtheme = theme_minimal(),palette = "npg") +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.position=c(0.13,0.78),
        legend.box.background = element_rect(color="grey", size=0.5),
       legend.text = element_text(size = 8))+
  labs(x = paste("PCoA 1 (", format(100*eig[1]/sum(eig), digits = 4), "%)",sep = ""), 
       y = paste("PCoA 2 (", format(100*eig[2]/sum(eig), digits = 4), "%)",sep = ""))

##legend function
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 7, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_blank(), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

tiff(filename = "Farm.AMR.ARO.PCoA.ggpubr.tif",width = 2500,height = 2000,res=600,compression="lzw")
addSmallLegend(p)
dev.off()


##########Procrustes analyses###################
#ARGs profile
AMR.gene <- read.table("adult425.faeces.ARO.xls",header=T,sep = "\t",row.names = 1,check.names=F)
AMR.gene$Farm <- factor(AMR.gene$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))
AMR.gene <- AMR.gene[order(AMR.gene$Farm),]
AMR.gene.Abun <- AMR.gene[,-c(1:5)]

#species data
species <- read.table("adult425.faeces.species.fpkm.xls",header = T,row.names=1,check.names = F,sep = "\t")

library(vegan)
#Hellinger transform
AMR.gene.Abun.hel <- decostand(AMR.gene.Abun, "hellinger")
species.data.hel <- decostand(species.data, "hellinger")

gene.distance <- vegdist(AMR.gene.Abun.hel,method = "bray") 
species.distance <- vegdist(species.data.hel,method = "bray")

# make pcoas 
library(ape)
pcoa_gene <- as.data.frame(pcoa(gene.distance)$vectors)
pcoa_species <- as.data.frame(pcoa(species.distance)$vectors)

# procrustes
pro <- procrustes(pcoa_gene, pcoa_species)
pro_test <- protest(pcoa_gene, pcoa_species, perm = 999)  
pro_test
# Call:
#   protest(X = pcoa_gene, Y = pcoa_species, permutations = 999) 
# 
# Procrustes Sum of Squares (m12 squared):        0.4293 
# Correlation in a symmetric Procrustes rotation: 0.7554 
# Significance:  0.001 
# 
# Permutation: free
# Number of permutations: 999

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100
beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "AMR gene"
beta_pro$Farm <- AMR.gene$Farm
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Species level"
trans_pro$Farm <- AMR.gene$Farm

colnames(trans_pro) <- colnames(beta_pro)

r <- signif(pro_test$t0, 3)
pval <- signif(pro_test$signif, 1)

plot <- rbind(beta_pro, trans_pro)

library(ggplot2)
library(ggsci)
plot$Farm <- factor(plot$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))
gene_species <- ggplot(plot) +
  geom_point(size = 0.5, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = Farm,shape = type)) +
  scale_color_npg() +   #point:scale_color_npg()£¬bar:scale_fill_npg()
  theme_bw()  +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName,color=Farm), alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'right',
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  xlab(paste0("PCoA 1 (",percent_var[1],"%)")) +
  ylab(paste0("PCoA 2 (",percent_var[2],"%)")) +
  labs(title=paste0("r = ",r,", p = ",pval),size=6)

addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 7, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize,face = "bold",vjust = 0), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
tiff(filename = "AMR.species.protest.bray.PCoA.tif",width = 2500,height = 2000,res=600,compression="lzw")
addSmallLegend(gene_species)
dev.off()

#The default color of ggplot2
library(scales)
#show_col(hue_pal()(9)) 
#cols <- c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
#cols <- hue_pal()(9)
cols <-pal_npg("nrc")(9)
show_col(cols)
table(beta_pro$Farm)
# Wild   KD-3800    KD-3480    KD-1400 NC-Tibetan      NC-F6    Dingnan 
# 6          7          6          8          6        293         63 
# Jiangying      Shahu 
# 16         20 
cols.data <- rep(cols,c(6,7,6,8,6,293,63,16,20))

##residuals plot
residuals <- residuals(pro_test)
tiff(filename = "gene.species.protest.bray.residuals.tif",width = 1800,height = 1500,res=300,compression="lzw")
plot(pro_test, kind=2,col=cols.data)#kind = 2 plots an impulse diagram of residuals.
legend("topright",inset = 0.02,box.col = "grey",
       legend = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"),
       col = cols,lty = 1,cex = 0.6)
dev.off()


######### at drug level ###############
Drug_num <- read.table("ARO.drug.num.xls",header = T,row.names = 1,check.names = F,sep = "\t")
Drug_num$Drug <- rownames(Drug_num)

Mechanism_num <- read.table("ARO.Mechanism_num.xls",header = T,check.names = F,sep = "\t")
library(reshape2)
Drug.melt <- melt(Drug_num,id.vars = "Drug")
Mechanism.melt <- melt(Mechanism_num,id.vars = "Mechanism")

library(ggplot2)
ggplot(Drug.melt,aes(x=variable,y=value,fill=Drug))+
  geom_bar(stat="identity") 
ggplot(Mechanism.melt,aes(x=variable,y=value,fill=Mechanism))+
  geom_bar(stat="identity") 
