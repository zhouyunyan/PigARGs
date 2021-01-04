######Indicators for the abundance and richness of resistance genes###########
ARO <- read.table("adult425.faeces.ARO.xls",header=T,sep="\t",check.names=F)

AMR_level <- read.delim("clipboard",header=T,sep="\t",check.names=F) #ARO.total.xlsx£ºARO
merge.data <- merge(ARO,AMR_level[,c(1,6,7)],by ="ID")

data <- merge.data[,7:355] #abundance of 349 ARGs

gene_name <- read.table("ARO_Freq_abun.xls",header=T,sep="\t",check.names=F)
data <- data[,match(gene_name$AROID,names(data))]
names(data) <- gene_name$ARO_name

Richness <- merge.data$ARO_Freq #total ARGs number of each number
Abun <- rowSums(data)  #total abundance of each sample


library(psych)
#abundance
single_Corr_Abun <- corr.test(data,Abun,use = "complete",method = "spearman")
single_Abun_r <- single_Corr_Abun$r
single_Abun_p <- single_Corr_Abun$p

#Richness
single_Corr_Richness <- corr.test(data,Richness,use = "complete",method = "spearman")
single_Richness_r <- single_Corr_Richness$r
single_Richness_p <- single_Corr_Richness$p

#the gene with the highest correction to the total AMR level was selected as initial gene subset.
single_mean_r <- (single_Abun_r+single_Richness_r)/2
single_mean_best_r <- max(single_mean_r)
#[1] 0.8522104
single_best_num <- which(single_mean_r== single_mean_best_r)
single_best_name <- rownames(single_mean_r)[single_best_num]

single_Abun_best_r <- single_Abun_r[single_best_num]
#[1] 0.8891706
single_Abun_best_p <- single_Abun_p[single_best_num]
#[1] 3.614059e-143
single_Rich_best_r <- single_Richness_r[single_best_num]
#[1] 0.8152503
single_Rich_best_p <- single_Richness_p[single_best_num]
#[1] 7.36329e-100

#the initial gene subset combined with every other gene and calculated the correction between the sum abundance of each pair gene and total AMR level to find the best combination of pair, and the gene subset was updated
library(psych)
select <- single_best_num
final_result <- c(single_best_name,single_Abun_max_r,single_Abun_best_p,single_Rich_max_r,single_Rich_best_p,single_mean_best_r)
for (j in 1:9) {
  result <- NULL
  seq <- seq(1,ncol(data),1)
  for (i in seq[-select]) {
    sum <- rowSums(data[c(select,i)])
    Abun_corr <- corr.test(sum,Abun,use = "complete",method = "spearman")
    Abun_cor <- Abun_corr$r
    Abun_p <- Abun_corr$p
    Rich_corr <- corr.test(sum,Richness,use = "complete",method = "spearman")
    Rich_cor <- Rich_corr$r
    Rich_p <- Rich_corr$p
    merge_r <- (Abun_cor+ Rich_cor)/2
    a <- c(names(data)[i],Abun_cor,Abun_p,Rich_cor,Rich_p,merge_r)
    result <- rbind(result,a)
  }
  result <- as.data.frame(result,stringsAsFactors = FALSE)
  names(result) <- c("gene","Abun_r","Abun_p","Richness_r","Richness_p","merge_r")

  result[,2:6] <- lapply(result[,2:6], as.numeric)
  r.best <- max(result$merge_r)
  r.name <- result$gene[which(result$merge_r== r.best)] 
  num <- which(names(data)== r.name)
  best.data <- result[which(result$merge_r== r.best),]
  final_result <- rbind(final_result,best.data)
  select <- c(select,num)
}
write.table(final_result,file = "predict_Abun_Richness.xls",row.names=F,sep = "\t")

#A line chart shows the best combination
final_result <- read.delim(file = "clipboard",header = T,sep = "\t",check.names = F) #abundance_richness.xlsx:abun_richness
final_result$r <-round(final_result$r ,digits = 3)
final_result$group <- factor(final_result$group, levels = c("Abundance","Richness","Combination"))
final_result$gene <- factor(final_result$gene, levels = final_result$gene[1:10])

final_result$ID <- as.factor(final_result$ID)

library(ggplot2)
tiff(filename = "abundance_richness.tif",width = 2800,height = 2500,res=600,compression="lzw")
ggplot(final_result, aes(x=final_result$gene, y=final_result$r, group=group,color=group)) + 
  geom_line(linetype="dotted") +
  geom_point(size=2, shape=20) +
  #geom_text(aes(label=gene),nudge_x = 0.15, nudge_y = -0.01,size = 2.8)+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.16),
        legend.box.background = element_rect(color="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.6))+
  labs(x="AMR genes",y="Spearman correlation")
dev.off()

final_result <- read.delim(file = "clipboard",header = T,sep = "\t",check.names = F) #abundance_richness.xlsx:predictor_best
final_result$r <-round(final_result$r ,digits = 3)
final_result$group <- factor(final_result$group, levels = c("Abundance","Richness","Combination"))
final_result$ID <- as.factor(final_result$ID)

tiff(filename = "predictor_best.tif",width = 3000,height = 2600,res=600,compression="lzw")
ggplot(final_result, aes(x=final_result$ID, y=final_result$r, group=group,color=group)) + 
  geom_line(linetype="dotted") +
  geom_point(size=2, shape=20) +
  geom_text(aes(label=gene),nudge_x = 0.15, nudge_y = -0.01,size = 2.6)+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.14),
        legend.box.background = element_rect(color="grey", size=0.5))+
   labs(x="Number of AMR genes",y="Spearman correlation")
dev.off()


###The prediction power of selected indicators for abundance and richness of resistance genes in other countries
data <- read.table(file = "stat.xls",header = T,row.names=1, sep = "\t",check.names = F)
data <- as.data.frame(t(data))

gene_Abun <- rowSums(data) #total abundance
gene_Abun <- as.data.frame(gene_Abun)

library(vegan)
gene_num <- specnumber(data) #total richness
gene_num <- as.data.frame(gene_num)


#correction for Richness
pig_data <- read.delim(file = "clipboard",header = T,sep = "\t",check.names = F) #abundance_richness.xlsx:Europe_Richness_select_8
sub_data <- pig_data[3:10]
rownames(sub_data) <- pig_data$ID

gene_Abun <- gene_Abun[rownames(sub_data),]
gene_num <- gene_num[rownames(sub_data),]

top10 <- read.table(file = "predict_richness_top10.xls",header = T,sep = "\t",check.names = F) 
names(sub_data) <- as.character(top10$gene[c(3:10)])

library(psych)
single_Corr_Rich <- corr.test(sub_data,gene_num,use = "complete",method = "spearman")
single_Rich_r <- single_Corr_Rich$r
# [,1]
# QnrS1        0.2723036
# aadA4       -0.1263220
# oqxB         0.5004521
# Erm(36)      0.2582825
# AAC(6')-IIa  0.1128578
# lsaB         0.4909802
# tet(B)       0.5109189
# lnuD         0.5383940
single_Rich_p <- single_Corr_Rich$p
# [,1]
# QnrS1       7.080869e-04
# aadA4       1.732819e-01
# oqxB        2.433755e-12
# Erm(36)     1.157311e-03
# AAC(6')-IIa 1.732819e-01
# lsaB        6.453915e-12
# tet(B)      7.574862e-13
# lnuD        2.158030e-14

#the cumulative correlation of gene subset:from sub_data[,1:2]to sub_data[,1:8]
pig_sum <- rowSums(sub_data[1:8])
pig_corr_Abun <- corr.test(pig_sum,gene_num,use = "complete",method = "spearman")
pig_corr_Abun$r
#0.6989048
pig_corr_Abun$p
#[1] 1.928173e-28

pig_data <- read.delim(file = "clipboard",header = T,sep = "\t",check.names = F) #abundance_richness.xlsx: Europe_combine_select6
sub_data <- sub_data[3:8]
rownames(sub_data) <- pig_data$ID

gene_Abun <- gene_Abun[match(rownames(sub_data),rownames(gene_Abun)),]
gene_num <- gene_num[match(rownames(sub_data),rownames(gene_num)),]

top10 <- read.table(file = "predict_Abun_Richness.xls",header = T,sep = "\t",check.names = F) 
names(sub_data) <- as.character(top10$gene[c(1,2,6,7,8,10)])


library(psych)
#both abundance and richness into account
single_Corr_Abun <- corr.test(sub_data,gene_Abun,use = "complete",method = "spearman")
single_Abun_r <- single_Corr_Abun$r
# [,1]
# aadA   0.48063777
# tetM   0.57208029
# vgaE   0.26084960
# catS   0.07537615
# QnrS1  0.20308066
# dfrA20 0.05972754
single_Abun_p <- single_Corr_Abun$p
# [,1]
# aadA   2.194748e-11
# tetM   1.079451e-16
# vgaE   1.342234e-03
# catS   6.157139e-01
# QnrS1  1.669218e-02
# dfrA20 6.157139e-01
single_Rich_p <- single_Corr_Rich$p

#Sum of 6 genes
pig_sum <- rowSums(sub_data)
pig_corr_Abun <- corr.test(pig_sum,gene_Abun,use = "complete",method = "spearman")
pig_corr_Abun$r
#0.5543573
pig_corr_Abun$p
#[1] 2.693926e-16

single_Corr_Rich <- corr.test(sub_data,gene_num,use = "complete",method = "spearman")
single_Rich_r <- single_Corr_Rich$r
# [,1]
# aadA   0.5765410
# tetM   0.7406747
# vgaE   0.6210510
# catS   0.1622451
# QnrS1  0.2723036
# dfrA20 0.3855345
single_Rich_p <- single_Corr_Rich$p
# [,1]
# aadA   3.548354e-17
# tetM   1.172737e-32
# vgaE   2.036289e-20
# catS   2.735221e-02
# QnrS1  3.540435e-04
# dfrA20 1.795321e-07

#Sum of 6 genes
pig_sum <- rowSums(sub_data)
pig_corr_Rich <- corr.test(pig_sum,gene_num,use = "complete",method = "spearman")
pig_corr_Rich$r
#0.7807094
pig_corr_Rich$p
#[1] 3.226565e-39

#the cumulative correlation of gene subset:from sub_data[,1:2] to sub_data[,1:6]
pig_sum <- rowSums(sub_data[,1:2])
pig_corr <- corr.test(pig_sum,gene_num,use = "complete",method = "spearman")
pig_corr$r

pig_corr$p


#The cumulative prediction power of the best gene subset that combines both the abundance and richness of the ARGs  for abundance  and richness  respectively. 
#abundance_richness.xlsx:Europe_best
final_result <- read.delim(file = "clipboard",header = T,sep = "\t",check.names = F)
final_result$r <-round(final_result$r ,digits = 3)
final_result$group <- factor(final_result$group, levels = c("Abundance","Richness","Combination"))
final_result$ID <- as.factor(final_result$ID)

library(ggplot2)

library(ggrepel) 
tiff(filename = "Europe_best.tif",width = 3000,height = 2600,res=600,compression="lzw")
ggplot(final_result, aes(x=final_result$ID, y=final_result$r, group=group,color=group)) + 
  geom_line(linetype="dotted") +
  geom_point(size=2, shape=20) +
  geom_text_repel(aes(label=gene),size = 2.6)+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.14),
        legend.box.background = element_rect(color="grey", size=0.5))+
  labs(x="Number of AMR genes",y="Spearman correlation")
dev.off()

final_result <- read.delim(file = "clipboard",header = T,sep = "\t",check.names = F) #abundance_richness.xlsx:Europe_best
final_result$r <-round(final_result$r ,digits = 3)
final_result$group <- factor(final_result$group, levels = c("Abundance","Richness","Combination"))
final_result$gene <- factor(final_result$gene,levels = final_result$gene[1:6])

library(ggplot2)
#install.packages("ggrepel")
library(ggrepel) 
tiff(filename = "Europe_Abun_Richness.tif",width = 3000,height = 2600,res=600,compression="lzw")
ggplot(final_result, aes(x=final_result$gene, y=final_result$r, group=group,color=group)) + 
  geom_line(linetype="dotted") +
  geom_point(size=2, shape=20) +
  #geom_text_repel(aes(label=gene),size = 2.6)+
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.14),
        legend.box.background = element_rect(color="grey", size=0.5))+
  labs(x="AMR genes",y="Spearman correlation")
dev.off()