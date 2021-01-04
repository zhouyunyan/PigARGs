###The distribution frequency and average abundance of 349 ARGs in all samples.

data <- read.table("ARO_Freq_abun.xls",header = T,check.names = F,sep = "\t")

data$Abun.log<- log10(data$Ave_Abun)

mean_abun <- mean(data$Abun.log)
data$group <- rep(NA,nrow(data))

data$group[which(data$Freq >=404 & data$Abun.log >= mean_abun)] <- "HFHA"
data$group[which(data$Freq >=404 & data$Abun.log < mean_abun)] <- "HFLA"
data$group[which(data$Freq < 404 & data$Abun.log >= mean_abun)] <- "LFHA"
data$group[which(data$Freq < 404 & data$Abun.log <  mean_abun)] <- "LFLA"

library(ggplot2)

p <- ggplot(data,aes(x=data$Freq,y=data$Abun.log,color=group))+ 
  geom_point(alpha=0.5,size=1)+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        legend.title=element_blank(),
        legend.position = c(0.12,0.88),
        legend.box.background = element_rect(color="grey", size=0.5))+
  labs(x="Number of samples",y="AMR gene abundance (log10)")+
  geom_vline(xintercept = 404,colour = "red") +geom_hline(yintercept =mean_abun,colour = "red")


#legend function
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 7, spaceLegend = 0.5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_blank(), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

tiff(filename = "ARO.Abun.sampleNum.tif",width = 2000,height = 1800,res=600,compression="lzw")
addSmallLegend(p)
dev.off()