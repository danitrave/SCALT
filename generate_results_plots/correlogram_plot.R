library("ComplexHeatmap")
library(colorRamp2)
library(corrplot)

#### Heatmap function ####
singleSampleHeatmap <- function(pathScalt,pathOriginal,NAME){
  #### SCALT Annotation ####
  df1 <- read.table(pathScalt,sep="\t",header = T)
  colnames(df1) <- c("SCALT")
  df1$SCALT <- toupper(df1$SCALT)
  cathegsSCALT <- data.frame(table(df1$SCALT))
  colnames(cathegsSCALT) <- c("CT","Freq")
  refinedCathegsSCALT <- cathegsSCALT[cathegsSCALT$Freq>100,]
  refinedCathegsSCALT <- refinedCathegsSCALT[!refinedCathegsSCALT$CT == "UNCLASSIFIED", ]
  
  #### Annotation original ####
  df2 <- read.table(pathOriginal,sep="\t",header=T)
  colnames(df2) <- c("original")
  df2$original <- toupper(df2$original)
  cathegsOriginal <- data.frame(table(df2$original))
  colnames(cathegsOriginal) <- c("CT","Freq")
  refinedCathegsOriginal <- cathegsOriginal[cathegsOriginal$Freq>100,]
  
  bound <- cbind(df2,df1)
  
  counting <- data.frame(table(bound[bound$original %in% refinedCathegsOriginal$CT & bound$SCALT %in% refinedCathegsSCALT$CT,]))
  contingency <- data.frame(matrix(data=0.0,nrow=nrow(refinedCathegsOriginal),ncol = nrow(refinedCathegsSCALT)))
  row.names(contingency)<-refinedCathegsOriginal$CT
  colnames(contingency)<-refinedCathegsSCALT$CT
  
  for (i in counting$original){
    for (j in counting$SCALT){
      contingency[i,j]<-counting[counting$original==i & counting$SCALT==j,"Freq"]
    }
  }

  bar1 <- HeatmapAnnotation(" " = anno_barplot(data.frame(apply(contingency, 2, sum)),
                                               height=unit(4,"cm"),add_numbers=T,numbers_offset=unit(0.5,"mm"),
                                               numbers_rot=90,gp=gpar(border ="black",fill="#5aae61",lty="blank")))
  bar2 <- rowAnnotation(" " = anno_barplot(data.frame(apply(contingency, 1, sum)),
                                           width=unit(4,"cm"),add_numbers=T,numbers_offset=unit(0.5,"mm"),
                                           gp=gpar(border ="black",fill="#5aae61",lty="blank")))
  PLOT <- Heatmap(contingency/apply(contingency, 1, sum),cluster_columns = F,cluster_rows = TRUE,show_row_dend = FALSE,
          show_column_dend = FALSE,row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 6),column_names_side = c("bottom"),
          row_names_side = c("left"),top_annotation=bar1,right_annotation = bar2,column_names_rot = 45,
          heatmap_legend_param = list(title = ""),col = colorRamp2(c(0, 0.5, 1),c("#f7f7f7","#c2a5cf","#762a83")))

  return(PLOT)
}

#### Correlogram function ####
correlogramPlot <- function(pathScalt,pathOriginal,NAME,CEXLAB,CEXLEG,NUMB){
  #### SCALT Annotation ####
  df1 <- read.table(pathScalt,sep="\t",header = T)
  colnames(df1) <- c("SCALT")
  df1$SCALT <- toupper(df1$SCALT)
  cathegsSCALT <- data.frame(table(df1$SCALT))
  colnames(cathegsSCALT) <- c("CT","Freq")
  refinedCathegsSCALT <- cathegsSCALT[cathegsSCALT$Freq>NUMB,]
  refinedCathegsSCALT <- refinedCathegsSCALT[!refinedCathegsSCALT$CT == "UNCLASSIFIED", ]
  
  #### Annotation original ####
  df2 <- read.table(pathOriginal,sep="\t",header=T)
  colnames(df2) <- c("original")
  df2$original <- toupper(df2$original)
  cathegsOriginal <- data.frame(table(df2$original))
  colnames(cathegsOriginal) <- c("CT","Freq")
  refinedCathegsOriginal <- cathegsOriginal[cathegsOriginal$Freq>NUMB,]
  
  bound <- cbind(df2,df1)
  
  counting <- data.frame(table(bound[bound$original %in% refinedCathegsOriginal$CT & bound$SCALT %in% refinedCathegsSCALT$CT,]))
  contingency <- data.frame(matrix(data=0.0,nrow=nrow(refinedCathegsOriginal),ncol = nrow(refinedCathegsSCALT)))
  row.names(contingency)<-refinedCathegsOriginal$CT
  colnames(contingency)<-refinedCathegsSCALT$CT
  
  for (i in counting$original){
    for (j in counting$SCALT){
      contingency[i,j]<-counting[counting$original==i & counting$SCALT==j,"Freq"]
    }
  }
  
  ratio_contingency <- as.matrix(contingency/apply(contingency, 1, sum))
  ratio_contingency[is.nan(ratio_contingency)] <- 0.0
  PLOT <- corrplot(as.matrix(ratio_contingency), method = 'circle',col.lim = c(0,1),#"YlOrRd"
           tl.srt = 45,col=COL1("Blues",200),tl.col = "black",tl.cex = CEXLAB,cl.cex = CEXLEG) # colorful number
  return(PLOT)
}

#### Plot the heatmaps ####
# pdf("zheng68k_heatmap.pdf")
# HP <- singleSampleHeatmap("zheng68k/SCALT_classification.tsv","zheng68k/labels.tsv","zheng68k")
# HP
# dev.off()
# pdf("zheng_sorted_heatmap.pdf")
# HP <- singleSampleHeatmap("zheng_sorted/SCALT_classification.tsv","zheng_sorted/labels.tsv","zheng_sorted")
# HP
# dev.off()
# pdf("Baron_Human_heatmap.pdf")
# HP <- singleSampleHeatmap("Baron_Human/SCALT_classification.tsv","Baron_Human/labels.tsv","Baron_Human")
# HP
# dev.off()
# pdf("Muraro_heatmap.pdf")
# HP <- singleSampleHeatmap("Muraro/SCALT_classification.tsv","Muraro/labels.tsv","Muraro")
# HP
# dev.off()
# pdf("Segerstolpe_heatmap.pdf")
# HP <- singleSampleHeatmap("Segerstolpe/SCALT_classification.tsv","Segerstolpe/labels.tsv","Segerstolpe")
# HP
# dev.off()

#### Plot the correlograms ####
## Baron Human
png("correlogram_Baron_Human.png",width = 1000, height = 1000)
CP <- correlogramPlot("Baron_Human/SCALT_classification.tsv","Baron_Human/labels.tsv","Baron_Human",1.5,2.0,100)
CP
dev.off()
## Zheng68K
png("correlogram_zheng68k.png",width = 1500, height = 1000)
CP <- correlogramPlot("zheng68k/SCALT_classification.tsv","zheng68k/labels.tsv","zheng68k",1.5,2.0,100)
CP
dev.off()
## Zheng sorted
png("correlogram_zheng_sorted.png",width = 1500, height = 1000)
CP <- correlogramPlot("zheng_sorted/SCALT_classification.tsv","zheng_sorted/labels.tsv","zheng68k",1.5,2.0,100)
CP
dev.off()
## Segerstolpe 
png("correlogram_Segerstolpe.png",width = 1000, height = 1000)
CP <- correlogramPlot("Segerstolpe/SCALT_classification.tsv","Segerstolpe/labels.tsv","zheng68k",1.5,2.0,100)
CP
dev.off()
## Muraro
png("correlogram_Muraro.png",width = 1000, height = 1000)
CP <- correlogramPlot("Muraro/SCALT_classification.tsv","Muraro/labels.tsv","zheng68k",1.5,2.0,100)
CP
dev.off()
## Cell lines 10Xv2 5 cl
png("correlogram_10X5CL.png",width = 1000, height = 1000)
CP <- correlogramPlot("10x_5cl/SCALT_classification.tsv","10x_5cl/labels.tsv","zheng68k",1.5,2.0,100)
CP
dev.off()
## Cell lines CELseq2 5cl
png("correlogram_CelSeq2_5cl.png",width = 1000, height = 1000)
CP <- correlogramPlot("CelSeq2_5cl/SCALT_classification.tsv","CelSeq2_5cl/labels.tsv","zheng68k",1.5,2.0,0)
CP
dev.off()

### PBMCs samples ###
## 10Xv2pbmc1
png("correlogram_10Xv2pbmc1.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/10Xv2pbmc1/SCALT_classification.tsv","pbmc/10Xv2pbmc1/10Xv2pbmc1_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## 10Xv2pbmc2
png("correlogram_10Xv2pbmc2.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/10Xv2pbmc2/SCALT_classification.tsv","pbmc/10Xv2pbmc2/10Xv2pbmc2_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## 10Xv3pbmc1
png("correlogram_10Xv3pbmc1.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/10Xv3pbmc1/SCALT_classification.tsv","pbmc/10Xv3pbmc1/10Xv3pbmc1_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## CEL-Seqpbmc1
png("correlogram_CEL-Seqpbmc1.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/CEL-Seqpbmc1/SCALT_classification.tsv","pbmc/CEL-Seqpbmc1/CEL-Seqpbmc1_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## CEL-Seqpbmc2
png("correlogram_CEL-Seqpbmc2.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/CEL-Seqpbmc2/SCALT_classification.tsv","pbmc/CEL-Seqpbmc2/CEL-Seqpbmc2_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## Drop-Seqpbmc1
png("correlogram_Drop-Seqpbmc1.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/Drop-Seqpbmc1/SCALT_classification.tsv","pbmc/Drop-Seqpbmc1/Drop-Seqpbmc1_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## Drop-Seqpbmc2
png("correlogram_Drop-Seqpbmc2.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/Drop-Seqpbmc2/SCALT_classification.tsv","pbmc/Drop-Seqpbmc2/Drop-Seqpbmc2_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## Seq-Wellpbmc1
png("correlogram_Seq-Wellpbmc1.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/Seq-Wellpbmc1/SCALT_classification.tsv","pbmc/Seq-Wellpbmc1/Seq-Wellpbmc1_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## Seq-Wellpbmc2
png("correlogram_Seq-Wellpbmc2.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/Seq-Wellpbmc2/SCALT_classification.tsv","pbmc/Seq-Wellpbmc2/Seq-Wellpbmc2_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## Smart-Seq2pbmc1
png("correlogram_Smart-Seq2pbmc1.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/Smart-Seq2pbmc1/SCALT_classification.tsv","pbmc/Smart-Seq2pbmc1/Smart-Seq2pbmc1_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## Smart-Seq2pbmc2
png("correlogram_Smart-Seq2pbmc2.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/Smart-Seq2pbmc2/SCALT_classification.tsv","pbmc/Smart-Seq2pbmc2/Smart-Seq2pbmc2_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## inDroppbmc1
png("correlogram_inDroppbmc1.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/inDroppbmc1/SCALT_classification.tsv","pbmc/inDroppbmc1/inDroppbmc1_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()
## inDroppbmc2
png("correlogram_inDroppbmc2.png",width = 2000, height = 1200)
CP <- correlogramPlot("pbmc/inDroppbmc2/SCALT_classification.tsv","pbmc/inDroppbmc2/inDroppbmc2_labels.tsv","zheng68k",1.5,1.5,0)
CP
dev.off()