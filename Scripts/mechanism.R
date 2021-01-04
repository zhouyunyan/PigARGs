########The abundance of each resistance mechanism per farm#####
ARO.mechanism <- read.table("ARO.mechanism.xls",header=T,check.names = F,sep="\t")
SampleInfo <- read.table("sample.425.ARO.Abun.xls",header = T,check.names = F,sep="\t")
data <- merge(SampleInfo,ARO.mechanism,by="ID",all.x = TRUE)
data <- data[,-6]

library(plyr)
data.Farm <- ddply(data,"Farm",numcolwise(mean))

library(reshape2)
data.melt1 <- melt(data.Farm,id.vars = "Farm")

data.melt1$Farm <- factor(data.melt1$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))
library(ggplot2)
library(ggsci)
library(RColorBrewer)
colourCount = length(unique(data.melt1$variable))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
p1 <- ggplot(data.melt1,aes(x=Farm,y=value,fill=variable))+
  geom_bar(stat="identity")+
  labs(y="AMR abundance (FPKM)",x=NULL,fill="Resistance Mechanism") + 
  scale_fill_manual(values = getPalette(colourCount))+
  theme_bw()+
  #ylim(0,2000)+
  #guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.6))
#legend function
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 4.5, spaceLegend = 0.2) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
tiff(filename = "Farm.Resistance_Mechanism.tif",width = 2800,height = 1800,res=600,compression="lzw")
addSmallLegend(p1)
dev.off()
