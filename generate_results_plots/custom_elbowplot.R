#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("pheatmap")
library(RColorBrewer)
library(stringr) 
library("ComplexHeatmap")
library(ggplot2)

data = read.table("../data_for_elbow_plot.tsv",header = T,sep="\t",row.names = 1)

png("elbow_plot.png",width= 10,height= 6,units= "in",res=1200)
ggplot(data=data, aes(x=Significance, y=Accuracy)) +geom_line(color="#02818a")+geom_point(color="#3690c0")+
  scale_x_continuous(breaks = c(0.009902183,0.018315639,0.030197383,0.049787068,0.082084999,
                                0.135335283,0.22313016,0.367879441,0.60653066,1.0),
                     labels = c("9","8","7","6","5","4","3","2","1","0"))+
  ggtitle("") +xlab("Likelihood Difference") + ylab("Accuracy")+theme_classic()
dev.off()
