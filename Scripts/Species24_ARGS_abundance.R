########the abundance variation of ARGs in 24 species in nine pig farms #########
species.data <- read.table(file = "farm_species_abun.xls",header=T,check.names=F,sep="\t")
library(reshape2)
melt <- melt(species.data,id.vars="species")
names(melt) <- c("species","Farm","Abundance")
melt$Farm <- factor(melt$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))


#Bubble chart
library(ggplot2)
mytheme <- theme_minimal()+
  theme(
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title=element_text(hjust =0.5),
    axis.line.y=element_line(linetype=1,color='grey'),
    axis.line.x=element_line(linetype=1,color='grey'),
    axis.ticks = element_line(linetype=2,color='grey'),
    panel.grid=element_line(linetype=2,color='grey'),
    legend.background = element_rect(fill="gray90", size=0,color='white'),
    legend.text=element_text(face="bold",size=10),
    legend.title=element_text(face="bold",size=8),
    axis.text=element_text(face="bold",size=10,colour = "black")
  )


tiff(filename = "species24.Farm.AMR_Abun.tif",width = 3500,height = 4000,res=600,compression="lzw")
ggplot(data=melt, mapping=aes(x=Farm,y=species))+
  geom_point(stat = "identity",aes(size=Abundance),alpha=0.7,show.legend = TRUE)+
  labs(x=NULL,y=NULL)+
  scale_size(range = c(1, 8))+
  mytheme
dev.off()



