#####Difference of resistome in different gut locations#########

data <- read.delim("clipboard",header = T,row.names=1,sep="\t",check.names=F) #factor_comparison.xlsx:SampleType_ARO_abun
head(data)

data.t <- t(data[,3:ncol(data)])
head(data.t)
dim(data.t)

data.t.filter <- data.t[rowSums(data.t)>0,] #delete ARO absent in all samples
dim(data.t.filter)

data.t.filter <- as.data.frame(data.t.filter,stringsAsFactors = F)
data.t.filter$AROID <- rownames(data.t.filter)
dim(data.t.filter)
head(data.t.filter)

ARO <- read.delim("clipboard",header = T,check.names = F) #factor_comparison.xlsx:ARO
head(ARO)

#add drug information to data.t.filter
data.ARO <- merge(data.t.filter,ARO[,c(1,3,4)],by="AROID")
rownames(data.ARO) <- data.ARO$ARO_rename
head(data.ARO)
class(data.ARO$AROID)
dim(data.ARO)

filter_data <- data[,c(1,2,which(colnames(data) %in% data.ARO$AROID))] #Extract genes that exist in at least one sample from the original data
filter_data$ID <- rownames(filter_data)
dim(filter_data)
head(filter_data)
library(reshape2)
md <- melt(filter_data, id=c("ID","SampleType","Farm"))
dim(md)
head(md)

library(ggpubr)
compare_F6 <- compare_means(value~SampleType, data=md[which(md$Farm=="NC-F6"),],p.adjust.method = "fdr",group.by="variable")
compare_Wild <- compare_means(value~SampleType, data=md[which(md$Farm=="Wild"),],p.adjust.method = "fdr",group.by="variable")

write.table(compare_F6,file="F6_SampleType_ARO.Wilcoxon.xls",row.names=F,quote=F,sep="\t")
write.table(compare_Wild,file="Wild_SampleType_ARO.Wilcoxon.xls",row.names=F,quote=F,sep="\t")


######### The comparison of abundance and number of ARGs in different gut locations ###########
#import data
data.Abun.num <- read.delim("clipboard",header = T,row.names=1,sep="\t",check.names=F) #factor_comparison.xlsx:location_Abun_Number

head(data.Abun.num)
data.Abun.num$SampleType <-  factor(data.Abun.num$SampleType,levels=c("cecum","faeces"))

##legend function
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 7, spaceLegend = 0.2) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

library(ggpubr)

compare_means(Abundance~SampleType, data=data.Abun.num,group.by = "Farm",p.adjust.method = "fdr")
compare_means(Number~SampleType, data=data.Abun.num,group.by = "Farm",p.adjust.method = "fdr")

SampleType.comparisons <- list(c("cecum","faeces"))

#boxplot
p1 <- ggboxplot(data.Abun.num[1:20,], x="SampleType", y="Abundance",fill = "SampleType",palette = "jco",add = "point",facet.by = "Farm")+  
  guides(fill=FALSE) +
  labs(y="AMR abundance (FPKM)",x=NULL) + 
  #geom_text_repel(aes(label=ID),size = 2.6)+  
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size = 14,angle = 45,vjust = 0.6)) +
  #facet_wrap(~Farm,scales="free") +   
  stat_compare_means(comparisons = SampleType.comparisons,label = "p.signif")
#dev.off()


p2 <- ggboxplot(data.Abun.num[21:nrow(data.Abun.num),], x="SampleType", y="Abundance",fill = "SampleType",palette = "jco",add = "point",facet.by = "Farm")+  
  guides(fill=FALSE) +
  labs(y="AMR abundance (FPKM)",x=NULL) + 
  #geom_text_repel(aes(label=ID),size = 2.6)+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size = 14,angle = 45,vjust = 0.6)) +
  #facet_wrap(~Farm,scales="free") +   
  stat_compare_means(comparisons = SampleType.comparisons,label = "p.signif")

p3 <- ggboxplot(data.Abun.num[1:20,], x="SampleType", y="Number",fill = "SampleType",palette = "jco",add = "point",facet.by = "Farm")+  
  guides(fill=FALSE) +
  labs(y="AMR gene numbers",x=NULL) +
  #geom_text_repel(aes(label=ID),size = 2.6)+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size = 14,angle = 45,vjust = 0.6)) +
  #facet_wrap(~Farm,scales="free") +   
  stat_compare_means(comparisons = SampleType.comparisons,label = "p.signif")

p4 <- ggboxplot(data.Abun.num[21:nrow(data.Abun.num),], x="SampleType", y="Number",fill = "SampleType",palette = "jco",add = "point",facet.by = "Farm")+  
  guides(fill=FALSE) +
  labs(y="AMR gene numbers",x=NULL) +
  #geom_text_repel(aes(label=ID),size = 2.6)+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size = 14,angle = 45,vjust = 0.6)) +
  #facet_wrap(~Farm,scales="free") +   
  stat_compare_means(comparisons = SampleType.comparisons,label = "p.signif")


p.merge1 <- ggarrange(p2,p4,ncol=2,nrow=1,widths = c(2,2))
p.merge2 <- ggarrange(p1,p3,ncol=2,nrow=1,widths = c(2,2))
p.merge1
p.merge2


#save figure
tiff(filename = "SampleType.Wild.tif",width = 2000,height = 2500,res=600,compression="lzw")
addSmallLegend(p.merge1)
dev.off()

tiff(filename = "SampleType.F6.tif",width = 2000,height = 2500,res=600,compression="lzw")
addSmallLegend(p.merge2)
dev.off()


######### The comparison of abundance between faeces and cecum samples at drug level ###########
#calculate abundance group at drug level
library(plyr)
ARO.drug <- ddply(data.ARO,"Drug_rename",numcolwise(sum))
rownames(ARO.drug) <- ARO.drug$Drug_rename
ARO.drug <- ARO.drug[,-1]
head(ARO.drug)
dim(ARO.drug)


ARO.drug.t <- as.data.frame(t(ARO.drug),stringsAsFactors = F)
ARO.drug.data <- cbind(data[,c(1,2)],ARO.drug.t)
ARO.drug.data$ID <- rownames(ARO.drug.data)
head(ARO.drug.data)
#write.table(ARO.drug.data,file = "SampleType_ARO.drug.xls",quote = F,sep = "\t",row.names = T)


library(reshape2)
md <- melt(ARO.drug.data, id=c("ID","SampleType","Farm"))
dim(md)
head(md)

library(ggpubr)
compare_F6 <- compare_means(value~SampleType, data=md[which(md$Farm=="NC-F6"),],p.adjust.method = "fdr",group.by="variable")
compare_Wild <- compare_means(value~SampleType, data=md[which(md$Farm=="Wild"),],p.adjust.method = "fdr",group.by="variable")

write.table(compare_F6,file="F6_SampleType_Drug.Wilcoxon.xls",row.names=F,quote=F,sep="\t")
write.table(compare_Wild,file="Wild_SampleType_Drug.Wilcoxon.xls",row.names=F,quote=F,sep="\t")

#Significantly different resistance in wild boar
sig_Wild_compare <- compare_Wild[which(compare_Wild$p.adj<=0.05),]
dim(sig_Wild_compare)

#Significantly different resistance in F6
sig_F6_compare <- compare_F6[which(compare_F6$p.adj<=0.05),]
dim(sig_F6_compare)


#show significant result using heatmap
row_num <- which(rownames(ARO.drug) %in% sig_F6_compare$variable)
sig_data <- ARO.drug[row_num,1:20]
sig_data.log10 <- log10(sig_data+1)

SampleType <- as.data.frame(data[,1]) #group information
rownames(SampleType)=rownames(data)
names(SampleType)[1] <- "SampleType"
SampleType$SampleType <- factor(SampleType$SampleType,levels=c("cecum","faeces"))

library(pheatmap)
p.sig <- pheatmap(sig_data.log10,
                  cluster_rows =T,cluster_cols =F,legend = T,
                  show_rownames=T,show_colnames=T,
                  annotation_col = SampleType,
                  gaps_col=10,border_color=NA)
p.sig


tiff(filename = "F6_Drug_SampleType_pheatmap.tif",width = 3800,height = 4000,res=600,compression="lzw")
p.sig
dev.off()

