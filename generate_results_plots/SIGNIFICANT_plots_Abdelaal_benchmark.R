library(ggplot2)
library(reshape2)
library(viridis)
library(ggsci)

#### INTRA WITHOUT MULTIASSIGNED ####
Methods <- c("SCALT","SVM-rejection","scPred","SVM","singleCellNet","ACTINN","CaSTLe","scmapcell","LDA",
            "scmapcluster","RF","SingleR","LAmbDA","NMC","CHETAH","scVI","scID","Cell_BLAST","kNN",
            "SCINA","DigitalCellSorter","Garnett-CV","Garnett-pretrained","Moana","Garnett-DE",
            "SCINA-DE","DigitalCellSorter-DE")
Samples <- c("Baron Human","Muraro","Segerstolpe","10X_5cl","CelSeq2_5cl","Zheng Sorted","Zheng68K")

Sample <- rep(Samples,length(Methods))

Method = c(rep("SCALT",length(Samples)),rep("SVM-rejection",length(Samples)),rep("scPred",length(Samples)),rep("SVM",length(Samples)),
           rep("singleCellNet",length(Samples)),rep("ACTINN",length(Samples)),rep("CaSTLe",length(Samples)),
           rep("scmapcell",length(Samples)),rep("LDA",length(Samples)),rep("scmapcluster",length(Samples)),
           rep("RF",length(Samples)),rep("SingleR",length(Samples)),rep("LAmbDA",length(Samples)),rep("NMC",length(Samples)),
           rep("CHETAH",length(Samples)),rep("scVI",length(Samples)),rep("scID",length(Samples)),rep("Cell_BLAST",length(Samples)),
           rep("kNN",length(Samples)),rep("SCINA",length(Samples)),rep("DigitalCellSorter",length(Samples)),
           rep("Garnett-CV",length(Samples)),rep("Garnett-pretrained",length(Samples)),rep("Moana",length(Samples)),
           rep("Garnett-DE",length(Samples)),rep("SCINA-DE",length(Samples)),rep("DigitalCellSorter-DE",length(Samples)))

Median_F1_score <- c(0.99,0.99,0.99,1.0,1.0,0.95,0.94,0.99,0.98,1.0,1.0,1.0,0.99,0.92,0.98,0.98,1.0,1.0,1.0,0.96,NA,0.98,0.97,1.0,1.0,1.0,0.95,0.7,
                     0.96,0.97,0.99,1.0,1.0,0.88,0.74,0.96,0.97,1.0,1.0,1.0,0.88,0.74,0.94,0.96,0.98,1.0,0.99,0.84,0.79,
                     0.98,0.97,1.0,1.0,1.0,0.73,0.64,0.97,0.96,0.99,1.0,1.0,0.63,0.66,0.95,0.97,1.0,1.0,1.0,0.73,0.44,
                     0.94,0.96,0.98,1.0,1.0,0.81,0.66,0.97,0.95,0.97,1.0,1.0,0.66,0.32,0.80,0.95,0.96,1.0,1.0,NA,0.4,
                     0.91,0.84,0.93,0.92,0.90,0.71,0.56,0.94,0.96,0.97,1.0,1.0,0.65,0.11,0.56,0.97,0.99,1.0,1.0,0.97,0.64,
                     0.59,0.95,0.85,1.0,1.0,0.61,0.42,0.89,0.79,0.08,1.0,0.99,0.91,0.74,0.96,0.95,0.85,1.0,0.98,0.45,0.54,
                     NA,NA,NA,NA,NA,1.0,1.0,NA,NA,NA,NA,NA,0.99,0.78,NA,NA,NA,NA,NA,0.94,0.60,NA,NA,NA,NA,NA,0.98,0.54,
                     NA,NA,NA,NA,NA,0.93,0.50,NA,NA,NA,NA,NA,0.65,0.37,NA,NA,NA,NA,NA,0.38,0.47,NA,NA,NA,NA,NA,0.0,0.0)


df <- data.frame(Method,Sample,Median_F1_score) #reorder(Method,value) 
factor_column = factor(df$Method,rev(unique(df$Method)))
df$block = factor_column
data <- melt(df)

plot <- ggplot(data, aes(x = Sample, y = block, fill = value)) +
  geom_tile() +geom_text(aes(label = value), color = "Black", size = 3) + 
  scale_fill_viridis(option="mako",discrete = FALSE,name="Median F1 score",na.value="#f0f0f0")+
  theme(legend.position="bottom",legend.box="horizontal")+
  labs(title = "",x = "",y = "") 
plot <- plot +guides(col = guide_legend(label.position = "top"))+guides(fill = guide_legend(title = "Median F1 score"))
plot = plot +scale_fill_bs5("green",reverse = T)
ggsave("significant_intra_median_F1_Abdelaal_benchmark.png",plot = plot,dpi = 1000)

#### INTRA - unclassified cells without multiassigned ####
Methods <- c("SCALT","SVM-rejection","scPred","SVM","singleCellNet","ACTINN","CaSTLe","scmapcell","LDA",
             "scmapcluster","RF","SingleR","LAmbDA","NMC","CHETAH","scVI","scID","Cell_BLAST","kNN",
             "SCINA","DigitalCellSorter","Garnett-CV","Garnett-pretrained","Moana","Garnett-DE",
             "SCINA-DE","DigitalCellSorter-DE")
Samples <- c("Baron Human","Muraro","Segerstolpe","10X_5cl","CelSeq2_5cl","Zheng Sorted","Zheng68K")

Sample <- rep(Samples,length(Methods))

Method = c(rep("SCALT",length(Samples)),rep("SVM-rejection",length(Samples)),rep("scPred",length(Samples)),rep("SVM",length(Samples)),
           rep("singleCellNet",length(Samples)),rep("ACTINN",length(Samples)),rep("CaSTLe",length(Samples)),
           rep("scmapcell",length(Samples)),rep("LDA",length(Samples)),rep("scmapcluster",length(Samples)),
           rep("RF",length(Samples)),rep("SingleR",length(Samples)),rep("LAmbDA",length(Samples)),rep("NMC",length(Samples)),
           rep("CHETAH",length(Samples)),rep("scVI",length(Samples)),rep("scID",length(Samples)),rep("Cell_BLAST",length(Samples)),
           rep("kNN",length(Samples)),rep("SCINA",length(Samples)),rep("DigitalCellSorter",length(Samples)),
           rep("Garnett-CV",length(Samples)),rep("Garnett-pretrained",length(Samples)),rep("Moana",length(Samples)),
           rep("Garnett-DE",length(Samples)),rep("SCINA-DE",length(Samples)),rep("DigitalCellSorter-DE",length(Samples)))

unlabelled <- c(0.0,0.0,0.0,0.0,0.0,14.2,10.0,1.5,1.6,1.9,0.0,0.0,23.5,61.8,10.8,8.5,10,0.4,1.1,41,NA,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                4.2,3.8,6.4,0.0,0.0,58.2,70.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,7.9,1.1,3.6,0.0,0.2,7.2,20.2,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,NA,0.0,
                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.9,1.1,0.1,0.0,0.1,4.9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                8.3,17.3,32.1,24.2,9.8,50.4,12,3.2,19.6,23.1,0.1,68.1,9.4,29,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                NA,NA,NA,NA,NA,1.9,14.1,NA,NA,NA,NA,NA,0.0,0.0,NA,NA,NA,NA,NA,33.5,70,NA,NA,NA,NA,NA,14.8,55.2,
                NA,NA,NA,NA,NA,0.0,0.0,NA,NA,NA,NA,NA,11,50.9,NA,NA,NA,NA,NA,1.6,5.3,NA,NA,NA,NA,NA,0.0,0.0)


df <- data.frame(Method,Sample,unlabelled)
factor_column = factor(df$Method,rev(unique(df$Method)))
df$block = factor_column
data <- melt(df)

plot <- ggplot(data, aes(x = Sample, y = block, fill = value)) +
  geom_tile() +geom_text(aes(label = value), color = "Black", size = 3) + 
  scale_fill_viridis(option="cividis",discrete = FALSE,name="Unlabeled (%)",na.value="#f0f0f0")+
  theme(legend.position="bottom",legend.box="horizontal")+
  labs(title = "",x = "",y = "") 
plot <- plot +guides(col = guide_legend(label.position = "top"))+guides(fill = guide_legend(title = "Unlabeled (%)"))#+theme(legend.margin=margin(-10, 0, 0, 0))
plot = plot + scale_fill_bs5("indigo")
ggsave("significant_intra_unlabelled_Abdelaal_benchmark.png",plot = plot,dpi = 1000)

#### Inter benchmark - cell lines  without multiassigned ####
CT_methods <- c("SCALT","SVM","scVI","SingleR","singleCellNet","scPred","scmapcluster","scmapcell","scID","Cell_BLAST","ACTINN",
                "SVM-rejection","LDA","RF","LAmbDA","CHETAH","kNN9","CaSTLe","NMC")
CT_samples <- c("10X_5cl","CelSe2_5cl")

Sample <- rep(CT_samples,length(CT_methods))

Method <- c(rep("SCALT",length(CT_samples)),rep("SVM",length(CT_samples)),rep("scVI",length(CT_samples)),rep("SingleR",length(CT_samples)),
            rep("singleCellNet",length(CT_samples)),rep("scPred",length(CT_samples)),rep("scmapcluster",length(CT_samples)),
            rep("scmapcell",length(CT_samples)),rep("scID",length(CT_samples)),rep("Cell_BLAST",length(CT_samples)),
            rep("ACTINN",length(CT_samples)),rep("SVM-rejection",length(CT_samples)),rep("LDA",length(CT_samples)),
            rep("RF",length(CT_samples)),rep("LAmbDA",length(CT_samples)),rep("CHETAH",length(CT_samples)),rep("kNN9",length(CT_samples)),
            rep("CaSTLe",length(CT_samples)),rep("NMC",length(CT_samples)))

Median_F1_score <- c(1.0,1.0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.99,0.98,1,0.98,0.99,1,0.96,0.90,0.98,0.91,0.84)

df <- data.frame(Method,Sample,Median_F1_score)
factor_column = factor(df$Method,rev(unique(df$Method)))
df$block = factor_column
data <- melt(df)

plot <- ggplot(data, aes(x = Sample, y = block, fill = value)) +#coord_fixed()+
  geom_tile() +geom_text(aes(label = value), color = "Black", size = 3) + scale_fill_viridis(option="mako",discrete = FALSE,name="Median F1 score")+
  theme(legend.position="bottom",legend.box="horizontal",axis.text.x = element_text(angle = 45,hjust = 1))+coord_fixed()+
  labs(title = "",x = "",y = "")
plot <- plot +guides(col = guide_legend(label.position = "top"))+guides(fill = guide_legend(title = "Median F1 score"))#+theme(legend.margin=margin(-10, 0, 0, 0))
plot = plot + scale_fill_bs5("green",reverse = T)
plot
ggsave("significant_inter_cell_lines_median_F1_Abdelaal_benchmark.png",plot = plot,dpi = 1000)

### Unclassified cell lines - inter benchamark without multiassigned ####
CT_methods <- c("SVM","scVI","SingleR","singleCellNet","scPred","scmapcluster","scmapcell","scID","Cell_BLAST","ACTINN",
                "SVM-rejection","LDA","RF","LAmbDA","CHETAH","kNN9","CaSTLe","NMC","SCALT")
CT_samples <- c("10X_5cl","CelSe2_5cl")

Sample <- rep(CT_samples,length(CT_methods))

Method <- c(rep("SCALT",length(CT_samples)),rep("SVM",length(CT_samples)),rep("scVI",length(CT_samples)),rep("SingleR",length(CT_samples)),
            rep("singleCellNet",length(CT_samples)),rep("scPred",length(CT_samples)),rep("scmapcluster",length(CT_samples)),
            rep("scmapcell",length(CT_samples)),rep("scID",length(CT_samples)),rep("Cell_BLAST",length(CT_samples)),
            rep("ACTINN",length(CT_samples)),rep("SVM-rejection",length(CT_samples)),rep("LDA",length(CT_samples)),
            rep("RF",length(CT_samples)),rep("LAmbDA",length(CT_samples)),rep("CHETAH",length(CT_samples)),rep("kNN9",length(CT_samples)),
            rep("CaSTLe",length(CT_samples)),rep("NMC",length(CT_samples)))

unlabelled <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.2,8.1,40.9,13.1,0.0,5.3,0.0,0.0,19.9,20.3,2.1,3.6,0.0,0.0,0.0,0.0,0.0,0.0,
                     0.0,0.0,0.0,0.0,0.4,0.2,0.0,0.0,0.0,0.0,0.0,0.0)

df <- data.frame(Method,Sample,unlabelled)
factor_column = factor(df$Method,rev(unique(df$Method)))
df$block = factor_column
data <- melt(df)

plot <- ggplot(data, aes(x = Sample, y = block, fill = value)) +#coord_fixed()+
  geom_tile() +geom_text(aes(label = value), color = "Black", size = 3) + scale_fill_viridis(option="cividis",discrete = FALSE,name="Unlabeled (%)")+
  theme(legend.position="bottom",legend.box="horizontal",axis.text.x = element_text(angle = 45,hjust = 1))+coord_fixed()+
  labs(title = "",x = "",y = "")
plot <- plot +guides(col = guide_legend(label.position = "top"))+guides(fill = guide_legend(title = "Unlabeled (%)"))#+theme(legend.margin=margin(-10, 0, 0, 0))
plot = plot + scale_fill_bs5("indigo")
ggsave("significant_inter_cell_lines_unlabeled_Abdelaal_benchmark.png",plot = plot,dpi = 1000)

### Inter benchmark - Pancreas datasets without multiassigned ####
PA_methods <- c("SCALT","scVI","SVM","ACTINN","scmapcell","SingleR","singleCellNet","RF","scID","CaSTLe","scPred","LAmbDA","CHETAH",
                "SVM-rejection","LDA","scmapcluster","Cell_BLAST","kNN","NMC")

PA_samples <- c("Baron Human","Muraro","Segerstolpe")

Sample <- rep(PA_samples,length(PA_methods))

Method <-  c(rep("SCALT",length(PA_samples)),rep("scVI",length(PA_samples)),rep("SVM",length(PA_samples)),rep("ACTINN",length(PA_samples)),rep("scmapcell",length(PA_samples)),
             rep("SingleR",length(PA_samples)),rep("singleCellNet",length(PA_samples)),rep("RF",length(PA_samples)),
             rep("scID",length(PA_samples)),rep("CaSTLe",length(PA_samples)),rep("scPred",length(PA_samples)),
             rep("LAmbDA",length(PA_samples)),rep("CHETAH",length(PA_samples)),rep("SVM-rejection",length(PA_samples)),
             rep("LDA",length(PA_samples)),rep("scmapcluster",length(PA_samples)),rep("Cell_BLAST",length(PA_samples)),
             rep("kNN",length(PA_samples)),rep("NMC",length(PA_samples)))

Median_F1_score <- c(0.99,1.0,1.0,0.99,0.98,1.0,0.99,0.97,1.0,0.98,0.97,1.0,0.99,0.98,0.95,0.93,0.97,0.99,0.92,0.95,1.0,0.91,0.96,0.91,
                     0.86,0.97,0.98,0.68,0.95,0.98,0.99,0.71,0.97,0.54,0.87,0.96,0.85,0.96,0.49,0.26,0.97,0.99,0.75,0.73,0.75,
                     1.0,0.32,1.0,0.29,0.55,0.91,0.27,0.91,0.44,0,0.44,0.02)

df <- data.frame(Method,Sample,Median_F1_score)
factor_column = factor(df$Method,rev(unique(df$Method)))
df$block = factor_column
data <- melt(df)

plot <- ggplot(data, aes(x = Sample, y = block, fill = value)) +#coord_fixed()+
  geom_tile() +geom_text(aes(label = value), color = "Black", size = 3) + scale_fill_viridis(option="mako",discrete = FALSE,name="Median F1 score")+
  theme(legend.position="bottom",legend.box="horizontal",axis.text.x = element_text(angle = 45,hjust = 1))+coord_fixed()+
  labs(title = "",x = "",y = "")
plot <- plot +guides(col = guide_legend(label.position = "top"))+guides(fill = guide_legend(title = "Median F1 score"))#+theme(legend.margin=margin(-10, 0, 0, 0))
plot = plot + scale_fill_bs5("green",reverse = T)
plot
ggsave("significant_inter_pancreas_median_F1_Abdelaal_benchmarkPlot.png",plot = plot,dpi = 1000)

### Unclassified pancreas - inter benchamrk without multiassigned ####
PA_methods <- c("SCALT","scVI","SVM","ACTINN","scmapcell","SingleR","singleCellNet","RF","scID","CaSTLe","scPred","LAmbDA","CHETAH",
                "SVM-rejection","LDA","scmapcluster","Cell_BLAST","kNN","NMC")

PA_samples <- c("Baron Human","Muraro","Segerstolpe")

Sample <- rep(PA_samples,length(PA_methods))

Method <-  c(rep("SCALT",length(PA_samples)),rep("scVI",length(PA_samples)),rep("SVM",length(PA_samples)),rep("ACTINN",length(PA_samples)),rep("scmapcell",length(PA_samples)),
             rep("SingleR",length(PA_samples)),rep("singleCellNet",length(PA_samples)),rep("RF",length(PA_samples)),
             rep("scID",length(PA_samples)),rep("CaSTLe",length(PA_samples)),rep("scPred",length(PA_samples)),
             rep("LAmbDA",length(PA_samples)),rep("CHETAH",length(PA_samples)),rep("SVM-rejection",length(PA_samples)),
             rep("LDA",length(PA_samples)),rep("scmapcluster",length(PA_samples)),rep("Cell_BLAST",length(PA_samples)),
             rep("kNN",length(PA_samples)),rep("NMC",length(PA_samples)))

unlabelled <- c(61.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,12.9,8.9,16.4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,29.1,15.3,49.2,
                0.0,0.0,0.0,52,38.2,42.2,0.0,0.0,0.0,3,0.8,7.6,0.0,0.0,0.0,0.0,0.0,0.0,88.5,32.7,36.7,75.1,14,9,
                0.0,0.0,0.0,0.0,0.0,0.0)

df <- data.frame(Method,Sample,unlabelled)
factor_column = factor(df$Method,rev(unique(df$Method)))
df$block = factor_column
data <- melt(df)

plot <- ggplot(data, aes(x = Sample, y = block, fill = value)) +#coord_fixed()+
  geom_tile() +geom_text(aes(label = value), color = "Black", size = 3) + scale_fill_viridis(option="cividis",discrete = FALSE,name="Unlabeled (%)")+
  theme(legend.position="bottom",legend.box="horizontal",axis.text.x = element_text(angle = 45,hjust = 1))+coord_fixed()+
  labs(title = "",x = "",y = "")
plot <- plot +guides(col = guide_legend(label.position = "top"))+guides(fill = guide_legend(title = "Unlabeled (%)"))#+theme(legend.margin=margin(-10, 0, 0, 0))
plot = plot + scale_fill_bs5("indigo")
plot
ggsave("significant_inter_unlabelled_pancreas_Abdelaal_benchmark.png",plot = plot,dpi = 1000)