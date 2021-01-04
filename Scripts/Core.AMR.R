##############Core resistome#############
data <- read.table("adult425.faeces.ARO.xls",header = T,sep = "\t",check.names = F)

#Average abundance of each farm
library(plyr)
Farm.mean <- ddply(data,"Farm",numcolwise(mean))
Farm.mean <- Farm.mean[,-2]
write.table(Farm.mean,file = "SampleToFarm.ARO.xls",row.names = F,quote = F,sep = "\t")

#select Core ARGs
Core <- read.delim("clipboard",header = T,check.names = F) #Core_ARGs.xlsx:Core_ARGs
Farm <- c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu")
Core.data <- Farm.mean[,match(Core$AROID,colnames(Farm.mean))]
rownames(Core.data) <- Farm.mean$Farm
colnames(Core.data) <- Core$name
Core.data <- Core.data[Farm,]

library(pheatmap)
library(gplots)
library(ggplot2)
Core.data <- as.data.frame(t(Core.data))
Mechanism <- as.data.frame(Core[,c(3,5)])
rownames(Mechanism) <- Core$name
names(Mechanism) <- c("Drug","Mechanism")
tiff(filename = "Core.AMR.pheatmap.tif",width = 4500,height = 4500,res=600,compression="lzw")
pheatmap(Core.data,scale = "row",cluster_rows =F,cluster_cols =F,
         color = bluered(100),annotation_row  = Mechanism,
         fontsize=11,annotation_legend=T)
dev.off()

#Legend is not displayed
tiff(filename = "Core.AMR.pheatmap1.tif",width = 3600,height = 4500,res=600,compression="lzw")
pheatmap(Core.data,scale = "row",cluster_rows =F,cluster_cols =F,
         color = bluered(100),annotation_row  = Mechanism,
         fontsize=11,annotation_legend=F)
dev.off()

#The average percent of the abundance of core ARGs to the total abundance of all ARGs per farm. 
data <- read.table("adult425.faeces.ARO.xls",header = T,sep = "\t",check.names = F)
Core <- read.delim("clipboard",header = T,check.names = F) #Core_ARGs.xlsx:Core_ARGs
Core.data <- data[,which(colnames(data) %in% Core$AROID)]

Core.sum <- apply(Core.data, 1, sum)
total.sum <- apply(data[,7:356], 1, sum)
Percent <- Core.sum/total.sum
result <- data.frame(data$ID,Core.sum,total.sum,Percent)
write.table(result,file = "Core.AMR.Percent.xls",row.names = F,quote = F,sep = "\t")

result$Farm <- data$Farm
Farm.mean <- tapply(result$Percent, result$Farm, mean)
Farm.nonCore <- 1-Farm.mean
Farm.data <- data.frame(names(Farm.mean),Farm.mean,Farm.nonCore)
names(Farm.data) <- c("Farm","Core","Non-core")
write.table(Farm.data,file = "Core.AMR.Farm.Percent.xls",row.names = F,quote = F,sep = "\t")

Farm.data <- read.table("Core.AMR.Farm.Percent.xls",header = T,sep = "\t",check.names = F)
library(reshape2)
Farm.data.melt <- melt(Farm.data,id.vars = "Farm")

#bar plot
Farm.data.melt$Farm <- factor(Farm.data.melt$Farm,levels = c("Wild","KD-3800","KD-3480","KD-1400","NC-Tibetan","NC-F6","Dingnan","Jiangying","Shahu"))

library(ggplot2)
library(ggsci)
library(RColorBrewer)
Farm.data.melt$variable <- factor(Farm.data.melt$variable,levels = c("Non-core","Core"))
p1 <- ggplot(Farm.data.melt,aes(x=Farm,y=value,fill=variable))+
  geom_bar(stat="identity")+
  labs(y="The percent of core AMR genes (%)",x=NULL,fill=NULL) + 
  theme_bw()+
  #ylim(0,2000)+
  #guides(fill=FALSE)+
  guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.6))

##legend function
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 7, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
tiff(filename = "Core.AMR.percent.tif",width = 2500,height = 2000,res=600,compression="lzw")
addSmallLegend(p1)
dev.off()
