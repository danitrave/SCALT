#!/usr/bin/env Rscript
library(RColorBrewer)
library(stringr) 
library(ggplot2)
library(viridis)
library(forcats)
library("ComplexHeatmap")
library("pheatmap")

##### INTRA without multiassigned ####
liu = read.csv("fu_benchmark_accuracies/Liu_clustering_comparison_5fold.csv",header = T,row.names = 1)
zheng68k = read.csv("fu_benchmark_accuracies/Zheng68K_clustering_comparison_5fold.csv",header = T,row.names = 1)
zheng_sorted = read.csv("fu_benchmark_accuracies/Zhengsort_clustering_comparison_5fold.csv",header = T,row.names = 1)
liu$SCALT = c(0.88,0.885,0.887,0.879,0.884)
zheng68k$SCALT = c(0.60,0.602,0.592,0.598,0.591)
zheng_sorted$SCALT = c(0.687,0.679,0.68,0.672,0.668)

liu_data = data.frame(sample=rep("Liu",ncol(liu)),
                      Accuracy=c(liu$Cell_blast,liu$ItClust,liu$TOSICA,liu$CHETAH,liu$SingleR,
                                 liu$scmapcell,liu$scmapcluster,liu$SVM,liu$scBERT,liu$scDeepSort,
                                 liu$SCALT),Tool=c(rep("Cell_blast",5),rep("ItClust",5),
                                                   rep("TOSICA",5),rep("CHETAH",5),rep("SingleR",5),
                                                   rep("scmapcell",5),rep("scmapcluster",5),
                                                   rep("SVM",5),rep("scBERT",5),rep("scDeepSort",5),
                                                   rep("SCALT",5)))

zheng68k_data = data.frame(sample=rep("Zheng68K",ncol(zheng68k)),
                           Accuracy=c(zheng68k$Cell_blast,zheng68k$ItClust,zheng68k$TOSICA,zheng68k$CHETAH,zheng68k$SingleR,
                                      zheng68k$scmapcell,zheng68k$scmapcluster,zheng68k$SVM,zheng68k$scBERT,zheng68k$scDeepSort,
                                      zheng68k$SCALT),Tool=c(rep("Cell_blast",5),rep("ItClust",5),
                                                        rep("TOSICA",5),rep("CHETAH",5),rep("SingleR",5),
                                                        rep("scmapcell",5),rep("scmapcluster",5),
                                                        rep("SVM",5),rep("scBERT",5),rep("scDeepSort",5),
                                                        rep("SCALT",5)))

zheng_sorted_data = data.frame(sample=rep("Zhengsort",ncol(zheng_sorted)),
                               Accuracy=c(zheng_sorted$Cell_blast,zheng_sorted$ItClust,zheng_sorted$TOSICA,zheng_sorted$CHETAH,zheng_sorted$SingleR,
                                          zheng_sorted$scmapcell,zheng_sorted$scmapcluster,zheng_sorted$SVM,zheng_sorted$scBERT,zheng_sorted$scDeepSort,
                                          zheng_sorted$SCALT),Tool=c(rep("Cell_blast",5),rep("ItClust",5),
                                                                 rep("TOSICA",5),rep("CHETAH",5),rep("SingleR",5),
                                                                 rep("scmapcell",5),rep("scmapcluster",5),
                                                                 rep("SVM",5),rep("scBERT",5),rep("scDeepSort",5),
                                                                 rep("SCALT",5)))

data = rbind(rbind(liu_data,zheng68k_data),zheng_sorted_data)

p<-ggplot(data, aes(x=reorder(sample,Accuracy,FUN=median,decreasing = T), y=Accuracy, fill=reorder(Tool,Accuracy,decreasing=T))) +
  geom_boxplot()+scale_fill_brewer(palette="Paired",name = "Method")+theme_classic()+
  theme(axis.title.x=element_blank())
png("intra_fu_without_multiassigned.png",width= 12,height= 8,units= "in",res=1200)
p
dev.off()

rm(list = ls())

##### INTRA with multiassigned #####
liu = read.csv("fu_benchmark_accuracies/Liu_clustering_comparison_5fold.csv",header = T,row.names = 1)
zheng68k = read.csv("fu_benchmark_accuracies/Zheng68K_clustering_comparison_5fold.csv",header = T,row.names = 1)
zheng_sorted = read.csv("fu_benchmark_accuracies/Zhengsort_clustering_comparison_5fold.csv",header = T,row.names = 1)
liu$SCALT = c(0.998,0.998,0.995,0.997,0.990)
zheng68k$SCALT = c(0.960,0.967,0.971,0.962,0.956)
zheng_sorted$SCALT = c(0.966,0.961,0.960,0.961,0.965)

liu_data = data.frame(sample=rep("Liu",ncol(liu)),
                      Accuracy=c(liu$Cell_blast,liu$ItClust,liu$TOSICA,liu$CHETAH,liu$SingleR,
                                 liu$scmapcell,liu$scmapcluster,liu$SVM,liu$scBERT,liu$scDeepSort,
                                 liu$SCALT),Tool=c(rep("Cell_blast",5),rep("ItClust",5),
                                                   rep("TOSICA",5),rep("CHETAH",5),rep("SingleR",5),
                                                   rep("scmapcell",5),rep("scmapcluster",5),
                                                   rep("SVM",5),rep("scBERT",5),rep("scDeepSort",5),
                                                   rep("SCALT",5)))

zheng68k_data = data.frame(sample=rep("Zheng68K",ncol(zheng68k)),
                           Accuracy=c(zheng68k$Cell_blast,zheng68k$ItClust,zheng68k$TOSICA,zheng68k$CHETAH,zheng68k$SingleR,
                                      zheng68k$scmapcell,zheng68k$scmapcluster,zheng68k$SVM,zheng68k$scBERT,zheng68k$scDeepSort,
                                      zheng68k$SCALT),Tool=c(rep("Cell_blast",5),rep("ItClust",5),
                                                             rep("TOSICA",5),rep("CHETAH",5),rep("SingleR",5),
                                                             rep("scmapcell",5),rep("scmapcluster",5),
                                                             rep("SVM",5),rep("scBERT",5),rep("scDeepSort",5),
                                                             rep("SCALT",5)))

zheng_sorted_data = data.frame(sample=rep("Zhengsort",ncol(zheng_sorted)),
                               Accuracy=c(zheng_sorted$Cell_blast,zheng_sorted$ItClust,zheng_sorted$TOSICA,zheng_sorted$CHETAH,zheng_sorted$SingleR,
                                          zheng_sorted$scmapcell,zheng_sorted$scmapcluster,zheng_sorted$SVM,zheng_sorted$scBERT,zheng_sorted$scDeepSort,
                                          zheng_sorted$SCALT),Tool=c(rep("Cell_blast",5),rep("ItClust",5),
                                                                     rep("TOSICA",5),rep("CHETAH",5),rep("SingleR",5),
                                                                     rep("scmapcell",5),rep("scmapcluster",5),
                                                                     rep("SVM",5),rep("scBERT",5),rep("scDeepSort",5),
                                                                     rep("SCALT",5)))

data = rbind(rbind(liu_data,zheng68k_data),zheng_sorted_data)

p<-ggplot(data, aes(x=reorder(sample,Accuracy,FUN=median,decreasing = T), y=Accuracy, fill=reorder(Tool,Accuracy,decreasing=T))) +
  geom_boxplot()+scale_fill_brewer(palette="Paired",name = "Method")+theme_classic()+
  theme(axis.title.x=element_blank())
png("intra_fu_with_multiassigned.png",width= 12,height= 8,units= "in",res=1200)
p
dev.off()

rm(list = ls())

#### INTER - scenario one ####
plot_heatmap <- function(d){
  PLOT1 = Heatmap(d,cluster_columns = F,cluster_rows = F,col = colorRampPalette(c("white","#e5d8bd","#cb181d"))(100),
                  show_column_names = T,show_row_names = T,show_heatmap_legend = F, column_names_rot=45,
                  na_col = "white",column_names_side="top",rect_gp = gpar(col = "grey60", lwd = 2),row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(d[i, j] > -1)
                      grid.text(sprintf("%.2f", d[i, j]), x, y, gp = gpar(fontsize = 18))
                  },top_annotation = HeatmapAnnotation(" " = anno_block(gp = gpar(fill = "#d9d9d9",col="#d9d9d9"))),
                  left_annotation = rowAnnotation("Testing dataset" = anno_block(gp = gpar(fill = "#d9d9d9",col="#d9d9d9"))),
                  row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 15))
  return(PLOT1)
}

unlabelled_heatmap <- function(d){
  PLOT1 = Heatmap(d,cluster_columns = F,cluster_rows = F,col = colorRampPalette(c("white","#c7e9b4","#41b6c4"))(100),
                  show_column_names = T,show_row_names = T,show_heatmap_legend = F, column_names_rot=45,
                  na_col = "white",column_names_side="top",rect_gp = gpar(col = "grey60", lwd = 2),row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(d[i, j] > -1)
                      grid.text(sprintf("%.2f", d[i, j]), x, y, gp = gpar(fontsize = 18))
                  },top_annotation = HeatmapAnnotation(" " = anno_block(gp = gpar(fill = "#d9d9d9",col="#d9d9d9"))),
                  left_annotation = rowAnnotation("Testing dataset" = anno_block(gp = gpar(fill = "#d9d9d9",col="#d9d9d9"))),
                  row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 15))
  return(PLOT1)
}

COLUMNS_NAMES = c("X10v2A_pbmc1","X10v2B_pbmc1","X10v3_pbmc1","DR_pbmc1","SW_pbmc1","iD_pbmc1","SM2_pbmc1","CL2_pbmc1")
ROW_NAMES = c("X10v2A_pbmc1","X10v2B_pbmc1","X10v3_pbmc1","DR_pbmc1","SW_pbmc1","iD_pbmc1","SM2_pbmc1","CL2_pbmc1")

## SVM ##
data = data.frame(c(0,0.88,0.92,0.83,0.80,0.82,0.68,0.64),c(0.91,0,0.93,0.82,0.68,0.58,0.50,0.70),
                  c(0.92,0.89,0,0.82,0.65,0.69,0.87,0.79),c(0.86,0.87,0.83,0,0.64,0.59,0.71,0.79),
                  c(0.90,0.81,0.86,0.74,0,0.83,0.79,0.84),c(0.92,0.92,0.93,0.87,0.75,0,0.92,0.88),
                  c(0.90,0.86,0.92,0.82,0.65,0.66,0,0.86),c(0.71,0.67,0.87,0.87,0.35,0.58,0.86,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_SVM.png",width= 14,height=10,units= "in",res=1200)
svm_p = plot_heatmap(data)
svm_p
dev.off()

## scDeepSort ##
data = data.frame(c(0,0.86,0.90,0.79,0.53,0.66,0.71,0.75),c(0.85,0,0.91,0.84,0.58,0.70,0.69,0.74),
                  c(0.79,0.84,0,0.79,0.74,0.56,0.82,0.73),c(0.41,0.52,0.65,0,0.33,0.47,0.60,0.71),
                  c(0.52,0.59,0.76,0.80,0,0.80,0.89,0.84),c(0.50,0.56,0.72,0.66,0.44,0,0.91,0.89),
                  c(0.66,0.70,0.80,0.69,0.54,0.79,0,0.77),c(0,0,0,0,0,0,0,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_scDeepSort.png",width= 14,height= 10,units= "in",res=1200)
scDeepSort_p = plot_heatmap(data)
scDeepSort_p
dev.off()

## scBERT ##
data = data.frame(c(0,0.87,0.92,0.77,0.75,0.59,0.88,0.85),c(0.88,0,0.93,0.74,0.55,0.31,0.90,0.84),
                  c(0.76,0.86,0,0.45,0.20,0.33,0.87,0.80),c(0.73,0.74,0.72,0,0.68,0.24,0.74,0.75),
                  c(0.87,0.74,0.82,0.75,0,0.78,0.93,0.77),c(0.39,0.33,0.34,0.40,0.36,0,0.47,0.38),
                  c(0.39,0.33,0.32,0.40,0.36,0.34,0,0.42),c(0,0,0,0,0,0,0,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_scBERT.png",width= 14,height= 10,units= "in",res=1200)
scBERT_p = plot_heatmap(data)
scBERT_p
dev.off()

## SingleR ##
data = data.frame(c(0,0.86,0.90,0.79,0.53,0.66,0.71,0.75),c(0.85,0,0.91,0.84,0.58,0.70,0.69,0.74),
                  c(0.79,0.84,0,0.79,0.74,0.56,0.82,0.73),c(0.41,0.52,0.65,0,0.33,0.47,0.60,0.71),
                  c(0.52,0.59,0.76,0.80,0,0.80,0.89,0.84),c(0.50,0.56,0.72,0.66,0.44,0,0.93,0.89),
                  c(0.66,0.70,0.80,0.69,0.54,0.79,0,0.77),c(0.57,0.61,0.64,0.64,0.67,0.58,0.77,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_SingleR.png",width= 14,height= 10,units= "in",res=1200)
SingleR_p = plot_heatmap(data)
SingleR_p
dev.off()

## ItClust ##
data = data.frame(c(0,0.71,0.75,0.62,0.57,0.66,0.53,0.57),c(0.83,0,0.81,0.72,0.51,0.59,0.71,0.45),
                  c(0.78,0.78,0,0.70,0.58,0.57,0.86,0.69),c(0.81,0.82,0.80,0,0.55,0.59,0.62,0.22),
                  c(0.90,0.89,0.91,0.81,0,0.78,0.89,0.87),c(0.84,0.82,0.83,0.75,0.66,0,0.86,0.76),
                  c(0.71,0.71,0.64,0.55,0.55,0.64,0,0.50),c(0.59,0.57,0.60,0.60,0.48,0.55,0.63,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_ItClust.png",width= 14,height= 10,units= "in",res=1200)
ItClust_p = plot_heatmap(data)
ItClust_p
dev.off()

## scmapcluster ##
data = data.frame(c(0,0.21,0.22,0.01,0.01,0,0,0),c(0.80,0,0.75,0.74,0.49,0.27,0.51,0.67),
                  c(0.79,0.80,0,0.67,0.11,0.41,0.73,0.66),c(0.52,0.65,0.71,0,0.12,0.02,0,0.09),
                  c(0.75,0.78,0.71,0.69,0,0.05,0,0.36),c(0,0,0.01,0,0,0,0,0),
                  c(0.03,0.02,0.13,0,0.01,0.01,0,0.03),c(0,0,0.03,0.03,0,0,0,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_scmapcluster.png",width= 14,height= 10,units= "in",res=1200)
scmapcluster_p = plot_heatmap(data)
scmapcluster_p
dev.off()

## TOSICA ##
data = data.frame(c(0,0.73,0.57,0.47,0.41,0.17,0.39,0.34),c(0.35,0,0.30,0.13,0.12,0.04,0.0,0.15),
                  c(0.62,0.68,0,0.53,0.12,0.43,0.34,0.26),c(0.45,0.56,0.63,0,0.35,0.22,0.70,0.66),
                  c(0.43,0.50,0.45,0.53,0,0.42,0.20,0.40),c(0.39,0.33,0.32,0.40,0.36,0,0.46,0.42),
                  c(0.07,0.07,0.09,0.01,0.01,0.03,0,0.07),c(0,0,0,0,0,0,0,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_TOSICA.png",width= 14,height= 10,units= "in",res=1200)
TOSICA_p = plot_heatmap(data)
TOSICA_p
dev.off()

## Cell BLAST ##
data = data.frame(c(0,0.74,0.71,0.62,0.60,0.59,0.16,0.26),c(0.64,0,0.70,0.47,0.47,0.46,0.23,0.36),
                  c(0.41,0.49,0,0.40,0.44,0.40,0.23,0.09),c(0.43,0.54,0.52,0,0.42,0.41,0.23,0.51),
                  c(0.59,0.45,0.48,0.39,0,0.48,0.12,0.17),c(0.56,0.45,0.53,0.42,0.49,0,0.15,0.11),
                  c(0.51,0.47,0.47,0.34,0.43,0.53,0,0.57),c(0.27,0.23,0.15,0.15,0.31,0.10,0.36,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_cell_blast.png",width= 14,height= 10,units= "in",res=1200)
cell_blast_p = plot_heatmap(data)
cell_blast_p
dev.off()

## CHETAH ##
data = data.frame(c(0,0.80,0.81,0.34,0.13,0.03,0.24,0.31),c(0.45,0,0.27,0.22,0.17,0.03,0.06,0.09),
                  c(0.55,0.68,0,0.24,0.05,0.16,0.16,0.17),c(0.59,0.67,0.78,0,0.30,0.17,0.34,0.57),
                  c(0.79,0.75,0.84,0.58,0,0.28,0.38,0.68),c(0.74,0.77,0.72,0.67,0.53,0,0.92,0.67),
                  c(0.43,0.49,0.71,0.40,0.33,0.33,0,0.59),c(0.43,0.48,0.56,0.55,0.30,0.26,0.52,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_CHETAH.png",width= 14,height= 10,units= "in",res=1200)
CHETAH_p = plot_heatmap(data)
CHETAH_p
dev.off()

## scmapcell ##
data = data.frame(c(0,0.44,0.33,0.39,0.24,0.07,0.26,0.24),c(0.65,0,0.63,0.60,0.48,0.50,0.45,0.40),
                  c(0.69,0.70,0,0.60,0.42,0.46,0.55,0.48),c(0.41,0.47,0.46,0,0.26,0.14,0.39,0.27),
                  c(0.59,0.52,0.55,0.50,0,0.47,0.54,0.46),c(0.23,0.14,0.15,0.08,0.21,0,0.18,0.12),
                  c(0.38,0.49,0.62,0.15,0.12,0.12,0,0.53),c(0.03,0.02,0.28,0.28,0.00,0.03,0.30,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenarion_1_scmapcell.png",width= 14,height= 10,units= "in",res=1200)
scmapcell_p = plot_heatmap(data)
scmapcell_p
dev.off()

## SCALT - without multiassigned ##
data = data.frame(c(0,0.82,0.84,0.86,0.76,0.85,0.71,0.70),c(0.87,0,0.88,0.79,0.72,0.90,0.74,0.78),
                  c(0.88,0.83,0,0.84,0.79,0.94,0.89,0.82),c(0.76,0.73,0.71,0,0.47,0.31,0.65,0.60),
                  c(0.93,0.85,0.88,0.87,0,0.97,0.80,0.84),c(0.90,0.79,0.80,0.94,0.87,0,0.87,0.68),
                  c(0.86,0.89,0.87,0.92,1,1,0,0.72),c(0.47,0.38,0.50,0.16,0.72,0.88,0.74,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenario_1_SCALT_without_multiassigned.png",width= 14,height= 10,units= "in",res=1200)
SCALT_p = plot_heatmap(data)
SCALT_p
dev.off()

## SCALT - with multiassigned ##
data = data.frame(c(0,0.94,0.94,0.96,0.90,0.99,0.73,0.84),c(0.98,0,0.96,0.96,0.92,1.0,0.78,0.91),
                  c(0.96,0.93,0,0.94,0.88,0.97,0.86,0.93),c(0.89,0.91,0.77,0,0.43,0.40,0.62,0.76),
                  c(0.99,0.97,0.98,0.99,0,0.99,0.86,0.99),c(0.97,0.93,0.94,0.98,0.93,0,0.97,0.98),
                  c(0.86,0.89,0.87,0.92,1,1,0,0.89),c(0.53,0.39,0.52,0.24,0.82,0.88,0.92,0))
colnames(data)=COLUMNS_NAMES
row.names(data)=ROW_NAMES
png("inter_scenario_1_SCALT_with_multiassigned.png",width= 14,height= 10,units= "in",res=1200)
SCALT_p = plot_heatmap(data)
SCALT_p
dev.off()

# ## SCALT - unlabeled without multiassigned ##
# unlabelled = data.frame(c(0,0.22,0.09,0.56,0.74,0.86,0.03,0.02),c(0.35,0,0.14,0.58,0.77,0.88,0.07,0.04),
#                         c(0.47,0.34,0,0.65,0.84,0.92,0.03,0.05),c(0.48,0.35,0.10,0,0.71,0.89,0.02,0.01),
#                         c(0.33,0.20,0.05,0.29,0,0.80,0.01,0.01),c(0.54,0.51,0.36,0.66,0.77,0,0.21,0.03),
#                         c(0.98,0.98,0.93,0.99,0.99,0.99,0,0.15),c(0.98,0.97,0.86,0.97,0.99,0.99,0.12,0))
# 
# colnames(unlabelled)=COLUMNS_NAMES
# row.names(unlabelled)=ROW_NAMES
# 
# ## SCALT - Unlabeled with multiassigned ##
# 
# unlabelled_SCALT_p = unlabelled_heatmap(unlabelled)
# png("inter_scenarion_1_unlabelled_SCALT.png",width= 14,height= 10,units= "in",res=1200)
# unlabelled_SCALT_p
# dev.off()

#### INTER - scenario two without multiassigned####
plot_heatmap_scenario_two <- function(d){
  PLOT1 = Heatmap(d,cluster_columns = F,cluster_rows = F,col = colorRampPalette(c("white","#e5d8bd","#cb181d"))(100),
                  show_column_names = T,show_row_names = T,show_heatmap_legend = F, column_names_rot=45,
                  na_col = "white",column_names_side="top",rect_gp = gpar(col = "grey60", lwd = 2),row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(d[i, j] > -1)
                      grid.text(sprintf("%.2f", d[i, j]), x, y, gp = gpar(fontsize = 15))
                  },row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 15))
  return(PLOT1)
}

scen2_COLUMNS_NAMES = c("SCALT","SVM","scDeepSort","scBERT","ItClust","scmapcluster","SingleR","CHETAH","TOSICA","Cell BLAST","scmapcell")
scen2_ROW_NAMES = c("X10v2A_pbmc1_pbmc2","X10v2B_pbmc1_pbmc2","DR_pbmc1_pbmc2","SW_pbmc1_pbmc2","iD_pbmc1_pbmc2","SM2_pbmc1_pbmc2","CL2_pbmc1_pbmc2")

data_scen2 = data.frame(c(0.87,0.83,0.85,0.90,0.86,0.86,0.58),c(0.91,0.86,0.86,0.79,0.88,0.88,0.41),c(0.83,0.84,0.71,0.72,0.68,0.84,0.0),
                        c(0.90,0.88,0.84,0.83,0.18,0.34,0.0),c(0.77,0.86,0.81,0.73,0.80,0.43,0.48),
                        c(0.09,0.86,0.51,0.30,0.0,0.15,0.0),c(0.83,0.84,0.71,0.72,0.68,0.84,0.49),
                        c(0.81,0.48,0.48,0.32,0.72,0.70,0.26),c(0.80,0.38,0.68,0.44,0.18,0.0,0.0),
                        c(0.73,0.56,0.76,0.40,0.49,0.64,0.34),c(0.48,0.72,0.61,0.49,0.32,0.40,0.0))

colnames(data_scen2)=scen2_COLUMNS_NAMES
row.names(data_scen2)=scen2_ROW_NAMES
scenario2_plot = plot_heatmap_scenario_two(data_scen2)
png("inter_scenarion_2_without_multiassigned.png",width= 12,height= 6,units= "in",res=1200)
scenario2_plot
dev.off()

unlabelled_scen2 = data.frame(c(0.24,0.27,0.16,0.69,0.76,0,0.01)) 

#### INTER - scenario two with multiassigned####
plot_heatmap_scenario_two <- function(d){
  PLOT1 = Heatmap(d,cluster_columns = F,cluster_rows = F,col = colorRampPalette(c("white","#e5d8bd","#cb181d"))(100),
                  show_column_names = T,show_row_names = T,show_heatmap_legend = F, column_names_rot=45,
                  na_col = "white",column_names_side="top",rect_gp = gpar(col = "grey60", lwd = 2),row_names_side = "left",
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    if(d[i, j] > -1)
                      grid.text(sprintf("%.2f", d[i, j]), x, y, gp = gpar(fontsize = 15))
                  },row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 15))
  return(PLOT1)
}

scen2_COLUMNS_NAMES = c("SCALT","SVM","scDeepSort","scBERT","ItClust","scmapcluster","SingleR","CHETAH","TOSICA","Cell BLAST","scmapcell")
scen2_ROW_NAMES = c("X10v2A_pbmc1_pbmc2","X10v2B_pbmc1_pbmc2","DR_pbmc1_pbmc2","SW_pbmc1_pbmc2","iD_pbmc1_pbmc2","SM2_pbmc1_pbmc2","CL2_pbmc1_pbmc2")

data_scen2 = data.frame(c(0.98,0.96,0.97,1.0,0.99,0.89,0.63),
                        c(0.91,0.86,0.86,0.79,0.88,0.88,0.41),c(0.83,0.84,0.71,0.72,0.68,0.84,0.0),
                        c(0.90,0.88,0.84,0.83,0.18,0.34,0.0),c(0.77,0.86,0.81,0.73,0.80,0.43,0.48),
                        c(0.09,0.86,0.51,0.30,0.0,0.15,0.0),c(0.83,0.84,0.71,0.72,0.68,0.84,0.49),
                        c(0.81,0.48,0.48,0.32,0.72,0.70,0.26),c(0.80,0.38,0.68,0.44,0.18,0.0,0.0),
                        c(0.73,0.56,0.76,0.40,0.49,0.64,0.34),c(0.48,0.72,0.61,0.49,0.32,0.40,0.0))

colnames(data_scen2)=scen2_COLUMNS_NAMES
row.names(data_scen2)=scen2_ROW_NAMES
scenario2_plot = plot_heatmap_scenario_two(data_scen2)
png("inter_scenarion_2_with_multiassigned.png",width= 12,height= 6,units= "in",res=1200)
scenario2_plot
dev.off()

unlabelled_scen2 = data.frame(c(0.24,0.27,0.16,0.69,0.76,0,0.01)) 