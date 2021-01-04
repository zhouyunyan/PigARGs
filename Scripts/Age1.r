#####Difference of resistome at different age#########
data <- read.delim("clipboard",header = T,row.names=1,sep="\t",check.names=F) factor_comparison.xlsx:Age_ARO_abun
head(data)

data.t <- t(data[,2:ncol(data)])
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
head(data.ARO)
class(data.ARO$AROID)
dim(data.ARO)

filter_data <- data[,c(1,which(colnames(data) %in% data.ARO$AROID))] #Extract genes that exist in at least one sample from the original data
filter_data$ID <- rownames(filter_data)
dim(filter_data)
head(filter_data)
library(reshape2)
md <- melt(filter_data, id=c("ID","Age"))
dim(md)
head(md)

library(ggpubr)
compare <- compare_means(value~Age, data=md,p.adjust.method = "fdr",group.by="variable")
head(compare)
write.table(compare,file="Age_ARO.Wilcoxon.xls",row.names=F,quote=F,sep="\t")

######### The comparison of abundance and number of ARGs in different gut locations ###########
data.Abun.num <- read.delim("clipboard",header = T,row.names=1,sep="\t",check.names=F) #factor_comparison.xlsx:Age_Abun_Number
head(data.Abun.num)
data.Abun.num$Age <-  factor(data.Abun.num$Age,levels=c("25d","240d"))

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
compare_means(Abundance~Age, data=data.Abun.num,p.adjust.method = "fdr")

Age.comparisons <- list(c("25d","240d"))

#abundance
p1 <- ggboxplot(data.Abun.num, x="Age", y="Abundance",fill = "Age",palette = "jco",add = "point")+  
  guides(fill=FALSE) +
  labs(y="AMR abundance (FPKM)",x=NULL) + 
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size = 14)) +
  stat_compare_means(comparisons = Age.comparisons,label = "p.signif")
p1

#Number
p2 <- ggboxplot(data.Abun.num, x="Age", y="Number",fill = "Age",palette = "jco",add = "point")+  
  guides(fill=FALSE) +
  labs(y="AMR gene numbers",x=NULL) + 
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size = 14)) +
  stat_compare_means(comparisons = Age.comparisons,label = "p.signif")
p2

#save
tiff(filename = "Age.drug.abun.com1.tif",width = 1000,height = 1600,res=600,compression="lzw")
addSmallLegend(p1)
dev.off()

tiff(filename = "Age.drug.number.com1.tif",width = 1000,height = 1600,res=600,compression="lzw")
addSmallLegend(p2)
dev.off()

######### The comparison of abundance between 25 and 240 days of age at drug level ###########
#calculate abundance group at drug level
library(plyr)
ARO.drug <- ddply(data.ARO,"Drug_rename",numcolwise(sum))
rownames(ARO.drug) <- ARO.drug$Drug_rename
ARO.drug <- ARO.drug[,-1]
head(ARO.drug)
dim(ARO.drug)


ARO.drug.t <- as.data.frame(t(ARO.drug),stringsAsFactors = F)
ARO.drug.data <- cbind(data$Age,ARO.drug.t)
names(ARO.drug.data)[1] <- "Age"
ARO.drug.data$ID <- rownames(ARO.drug.data)
head(ARO.drug.data)
#write.table(ARO.drug.data,file = "Age_ARO.drug.xls",quote = F,sep = "\t",row.names = T)

library(reshape2)
md <- melt(ARO.drug.data, id=c("ID","Age"))
dim(md)
head(md)

library(ggpubr)
compare_F6 <- compare_means(value~Age, data=md,p.adjust.method = "fdr",group.by="variable")
write.table(compare_F6,file="F6_Age_Drug.Wilcoxon.xls",row.names=F,quote=F,sep="\t")

#Significantly different resistance
sig_F6_compare <- compare_F6[which(compare_F6$p.adj<=0.05),]
dim(sig_F6_compare)

show significant result using heatmap
row_num <- which(rownames(ARO.drug) %in% sig_F6_compare$variable)
sig_data <- ARO.drug[row_num,1:20]
sig_data.log10 <- log10(sig_data+1)

Age <- as.data.frame(data[,1]) #group information
rownames(Age)=rownames(data)
names(Age)[1] <- "Age"
Age$Age <- factor(Age$Age,levels=c("25d","240d"))

library(pheatmap)
p.sig <- pheatmap(sig_data.log10,
               cluster_rows =T,cluster_cols =F,legend = T,
               show_rownames=T,show_colnames=T,
               annotation_col = Age,
               gaps_col=10,border_color=NA)
p.sig

tiff(filename = "F6_Drug_Age_pheatmap.tif",width = 3800,height = 4000,res=600,compression="lzw")
p.sig
dev.off()
