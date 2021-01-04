
species <- read.table("species.fpkm.3530.txt",header = T,row.names=1,check.names = F,sep = "\t")
ARO.data <- read.table("adult425.faeces.ARO.xls",header = T,sep = "\t",check.names = F)
sampleID <- ARO.data$ID
CARD_species <- read.table("clipboard",header = F,check.names = F)#CARD_species.xlsx:species_Freq_abundance

species_data <- species[rownames(species) %in% sampleID,colnames(species) %in% CARD_species$V1]
write.table(species_data,file = "CARD_species.Abun.xls",quote = F,sep = "\t")

species_data$Farm <- ARO.data$Farm

library(plyr)
species.Farm.mean <- ddply(species_data,"Farm",numcolwise(mean))
write.table(species.Farm.mean,file = "species.Farm.mean.xls",row.names = F,quote = F,sep = "\t")

species.Farm.mean <- read.table(file = "species.Farm.mean.xls",header = T,sep = "\t")
library(reshape2)
melt <- melt(species.Farm.mean,id.vars="Farm")
melt$Farm <- factor(melt$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))
melt$variable <- as.character(melt$variable)
melt <- melt[order(melt$variable),]
#Dingnan
melt$value[which(melt$value == max(melt$value))] <- 5000
names(melt)[3] <- "Abundance"

#Bubble graph
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

  
tiff(filename = "species24.Farm2.tif",width = 3500,height = 4000,res=600,compression="lzw")
#ggplot(data=melt, mapping=aes(x=Farm,y=variable,color=Abundance))+
ggplot(data=melt, mapping=aes(x=Farm,y=variable))+ 
  geom_point(stat = "identity",aes(size=Abundance),alpha=0.7,show.legend = TRUE)+
  labs(x=NULL,y=NULL)+
  scale_size(range = c(1, 8))+
  mytheme
#scale_color_manual(values=c("#666666","#FF0016"))
dev.off()