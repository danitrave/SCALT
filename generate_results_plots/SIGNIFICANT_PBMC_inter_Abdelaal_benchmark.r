## PBMC across dataset ----
#setwd('Benchmarking_paper/Classifiers/benchmarking')
source("Figure_scripts/evaluate.r")
library(tidyr)
library(tibble)
library(ggplot2)
library(viridis)
library(ggsci)

## Function to gather the data into key-value pairs, because that is what ggplot2 wants as input ----
gather_data2 <- function(data,classifier,train,test,exp){
  
  data <- data.frame(t(data))
  names(data) <- classifier
  score <- apply(data,2,mean)
  data$Train = train
  data$Test = test
  data$Exp = exp 
  data$Exp <- factor(data$Exp, 
                             levels = exp, ordered = TRUE)
  data = gather(data, "Classifier", "Value", -Train, -Test, -Exp, convert = TRUE)
  
  result <- list(data = data, score = score)
  
  return(result)
}

gather_data <- function(data,classifier,exp){
  
  data <- data.frame(data, row.names = classifier)
  names(data) <- exp
  data = rownames_to_column(data, var = "Classifier")
  data = gather(data, "Dataset", "Value", -Classifier, convert = TRUE)
  
  return(data)
}


## Read the data ----

# Indices per dataset
Ind_SM2 <- c(1,5890,8880,9123,12086,15298,18411,18634) #last index is to define the end of the dataset
Ind_10XV2 <- c(1,254,3476,3729,6951,10173,13349,16711)
Ind_10XV3 <- c(1,254,6660,6913,10107,13319,16469)
Ind_CL <- c(1,254,6572,9756,12904,16116,19229,19502)
Ind_DR <- c(1,254,6698,9920,10173,13395,16571,19933)
Ind_iD <- c(1,254,6181,9171,9414,12405,15544,18611)
Ind_SW <- c(1,246,6086,9016,9238,12126,15240,15791)

# Indices for Garnett CV
Gar_SM2 <- c(1,5457, 8177,8402,11338,14446,17521,17744) #last index is to define the end of the dataset
Gar_10XV2 <- c(1,240,3192,3427,6594,9702,12814,16114)
Gar_10XV3 <- c(1,240,6213,6448,9615,12723,15835)
Gar_CL <- c(1,240,6125,9039,12160,15268,18343,18616)
Gar_DR <- c(1,240,6213,9165,9400,12508,15620,18982)
Gar_iD <- c(1,240,5696,8416,8641,11577,14652,17646)
Gar_SW <- c(1,232,5601,8261,8465,11298,14298,14849)

# Lambda also needs special indices
Train_Idx <- c(253, 6444, 3222, 253, 3222, 3222, 3176)

# Name of the folder
Classifiers <- c('CaSTLe','CHETAH','scID','scmapcell','scmapcluster','scPred','singleCellNet','SingleR',
                 'ACTINN', 'Cell_BLAST', 'kNN9', 'LAmbDA', 'LDA', 'NMC', 'RF', 'scVI', 'SVM','SVM rejection','Garnett_CV',
                 'Garnett','SCINA','DigitalCellSorter')

# Name of the file
Filename <- c('CaSTLe','CHETAH','scID','scmapcell','scmapcluster','scPred','singleCellNet','SingleR',
               'ACTINN', 'Cell_BLAST', 'kNN9', 'LAmbDA', 'LDA', 'NMC', 'RF', 'scVI', 'SVM','SVM','Garnett_CV',
               'Garnett_5','SCINA_20','DigitalCellSorter_5')

Names <- c('CaSTLe','CHETAH','scID','scmapcell','scmapcluster','scPred','singleCellNet','SingleR',
           'ACTINN', 'Cell_BLAST', 'kNN9', 'LAmbDA', 'LDA', 'NMC', 'RF', 'scVI', 'SVM','SVM rejection','Garnett_CV',
           'Garnett_5','SCINA_20', 'DigitalCellSorter_5')

# How the name should appear in the figure
Expr <- c('CaSTLe','CHETAH','scID','scmapcell','scmapcluster','scPred','singleCellNet','SingleR',
          'ACTINN', 'Cell_BLAST', 'kNN9', 'LAmbDA', 'LDA', 'NMC', 'RF', 'scVI', 'SVM',expression('SVM'[rejection]),
          expression('Garnett'[CV]), expression('Garnett'['5']),expression('SCINA'['20']), expression('DigitalCellSorter'['5']))

# One heatmap
Exp_Big <- c('10XV2 - SM2', '10XV2 - 10XV3', '10XV2 - CL', '10XV2 - DR', '10XV2 - iD', '10XV2 - SW', '10XV2 - 10XV2',
             '10XV3 - SM2', '10XV3 - 10XV2', '10XV3 - CL', '10XV3 - DR', '10XV3 - iD', '10XV3 - SW',
             'CL - SM2', 'CL - 10XV2', 'CL - 10XV3', 'CL - DR', 'CL - iD', 'CL - SW', 'CL - CL',
             'DR - SM2', 'DR - 10XV2', 'DR - 10XV3', 'DR - CL', 'DR - iD', 'DR - SW', 'DR - DR',
             'iD - SM2', 'iD - 10XV2', 'iD - 10XV3', 'iD - CL', 'iD - DR', 'iD - SW', 'iD - iD',
             'SM2 - 10XV2', 'SM2 - 10XV3', 'SM2 - CL', 'SM2 - DR', 'SM2 - iD', 'SM2 - SW', 'SM2 - SM2',
             'SW - SM2', 'SW - 10XV2', 'SW - 10XV3', 'SW - CL', 'SW - DR', 'SW - iD', 'SW - SW')

Train <- c('10XV2', '10XV2', '10XV2', '10XV2', '10XV2', '10XV2', '10XV2',
             '10XV3', '10XV3', '10XV3', '10XV3', '10XV3', '10XV3',
             'CL', 'CL', 'CL', 'CL', 'CL', 'CL', 'CL',
             'DR', 'DR', 'DR', 'DR', 'DR', 'DR', 'DR',
             'iD', 'iD', 'iD', 'iD', 'iD', 'iD', 'iD',
             'SM2', 'SM2', 'SM2', 'SM2', 'SM2', 'SM2', 'SM2',
             'SW', 'SW', 'SW', 'SW', 'SW', 'SW', 'SW')

Test <- c('SM2', '10XV3', 'CL', 'DR', 'iD', 'SW', '10XV2',
             'SM2', '10XV2', 'CL', 'DR', 'iD', 'SW',
             'SM2', '10XV2', '10XV3', 'DR', 'iD', 'SW', 'CL',
             'SM2', '10XV2', '10XV3', 'CL', 'iD', 'SW', 'DR',
             'SM2', '10XV2', '10XV3', 'CL', 'DR', 'SW', 'iD',
          '10XV2', '10XV3', 'CL', 'DR', 'iD', 'SW', 'SM2',
             'SM2', '10XV2', '10XV3', 'CL', 'DR', 'iD', 'SW')


MedF1_Big <- matrix(data = NA, nrow = length(Classifiers), ncol = 48)
Unl_Big <- matrix(data = NA, nrow = length(Classifiers), ncol = 48)


for (j in c(1:(length(Classifiers)))){
  Count_Big <- 1
  Count_Small <- 1
  

  for (i in c(1:(length(Ind_10XV2)-1))){
    truefile <- paste("results/PBMC_benchmarking/10Xv2/",Classifiers[j],"/",Filename[j],"_True_Labels.csv",sep="")
    predfile <- paste("results/PBMC_benchmarking/10Xv2/",Classifiers[j],"/",Filename[j],"_Pred_Labels.csv",sep="")
    
    if (Classifiers[j] == 'LAmbDA'){
      res <- evaluate(truefile, predfile, Indices = c((Ind_10XV2[i]+Train_Idx[2]):(Ind_10XV2[i+1]+Train_Idx[2]-1)))
    } else if(Classifiers[j] == 'Garnett_CV'){
      res <- evaluate(truefile, predfile, Indices = c(Gar_10XV2[i]:(Gar_10XV2[i+1]-1)))
    } else {
      res <- evaluate(truefile, predfile, Indices = c(Ind_10XV2[i]:(Ind_10XV2[i+1]-1)))
    }
    
      MedF1_Big[j,Count_Big] <- res$MedF1
      Unl_Big[j,Count_Big] <- res$PercUnl
      Count_Big <- Count_Big + 1
  }
  
  for (i in c(1:(length(Ind_10XV3)-1))){
    truefile <- paste("results/PBMC_benchmarking/10Xv3/",Classifiers[j],"/",Filename[j],"_True_Labels.csv",sep="")
    predfile <- paste("results/PBMC_benchmarking/10Xv3/",Classifiers[j],"/",Filename[j],"_Pred_Labels.csv",sep="")
    
    if (Classifiers[j] == 'LAmbDA'){
      res <- evaluate(truefile, predfile, Indices = c((Ind_10XV3[i]+Train_Idx[3]):(Ind_10XV3[i+1]+Train_Idx[3]-1)))
    } else if(Classifiers[j] == 'Garnett_CV'){
      res <- evaluate(truefile, predfile, Indices = c(Gar_10XV3[i]:(Gar_10XV3[i+1]-1)))
    } else {
      res <- evaluate(truefile, predfile, Indices = c(Ind_10XV3[i]:(Ind_10XV3[i+1]-1)))
    }
    
    MedF1_Big[j,Count_Big] <- res$MedF1
    Unl_Big[j,Count_Big] <- res$PercUnl
    Count_Big <- Count_Big + 1
    
  }
  
  for (i in c(1:(length(Ind_CL)-1))){
    truefile <- paste("results/PBMC_benchmarking/CEL/",Classifiers[j],"/",Filename[j],"_True_Labels.csv",sep="")
    predfile <- paste("results/PBMC_benchmarking/CEL/",Classifiers[j],"/",Filename[j],"_Pred_Labels.csv",sep="")
    
    if (Classifiers[j] == 'LAmbDA'){
      res <- evaluate(truefile, predfile, Indices = c((Ind_CL[i]+Train_Idx[4]):(Ind_CL[i+1]+Train_Idx[4]-1)))
    } else if(Classifiers[j] == 'Garnett_CV'){
      res <- evaluate(truefile, predfile, Indices = c(Gar_CL[i]:(Gar_CL[i+1]-1)))
    } else {
      res <- evaluate(truefile, predfile, Indices = c(Ind_CL[i]:(Ind_CL[i+1]-1)))
    }
    
      MedF1_Big[j,Count_Big] <- res$MedF1
      Unl_Big[j,Count_Big] <- res$PercUnl
      Count_Big <- Count_Big + 1

  }
  
  for (i in c(1:(length(Ind_DR)-1))){
    truefile <- paste("results/PBMC_benchmarking/Drop/",Classifiers[j],"/",Filename[j],"_True_Labels.csv",sep="")
    predfile <- paste("results/PBMC_benchmarking/Drop/",Classifiers[j],"/",Filename[j],"_Pred_Labels.csv",sep="")
    
    if (Classifiers[j] == 'LAmbDA'){
      res <- evaluate(truefile, predfile, Indices = c((Ind_DR[i]+Train_Idx[5]):(Ind_DR[i+1]+Train_Idx[5]-1)))
    } else if(Classifiers[j] == 'Garnett_CV'){
      res <- evaluate(truefile, predfile, Indices = c(Gar_DR[i]:(Gar_DR[i+1]-1)))
    } else {
      res <- evaluate(truefile, predfile, Indices = c(Ind_DR[i]:(Ind_DR[i+1]-1)))
    }
    
      MedF1_Big[j,Count_Big] <- res$MedF1
      Unl_Big[j,Count_Big] <- res$PercUnl
      Count_Big <- Count_Big + 1

  }
  
  for (i in c(1:(length(Ind_iD)-1))){
    truefile <- paste("results/PBMC_benchmarking/inDrop/",Classifiers[j],"/",Filename[j],"_True_Labels.csv",sep="")
    predfile <- paste("results/PBMC_benchmarking/inDrop/",Classifiers[j],"/",Filename[j],"_Pred_Labels.csv",sep="")
    
    if (Classifiers[j] == 'LAmbDA'){
      res <- evaluate(truefile, predfile, Indices = c((Ind_iD[i]+Train_Idx[6]):(Ind_iD[i+1]+Train_Idx[6]-1)))
    } else if(Classifiers[j] == 'Garnett_CV'){
      res <- evaluate(truefile, predfile, Indices = c(Gar_iD[i]:Gar_iD[i+1]-1))
    } else {
      res <- evaluate(truefile, predfile, Indices = c(Ind_iD[i]:Ind_iD[i+1]-1))
    }
    
      MedF1_Big[j,Count_Big] <- res$MedF1
      Unl_Big[j,Count_Big] <- res$PercUnl
      Count_Big <- Count_Big + 1

  }
  
  for (i in c(1:(length(Ind_SM2)-1))){
    truefile <- paste("results/PBMC_benchmarking/SM2/",Classifiers[j],"/",Filename[j],"_True_Labels.csv",sep="")
    predfile <- paste("results/PBMC_benchmarking/SM2/",Classifiers[j],"/",Filename[j],"_Pred_Labels.csv",sep="")
    
    if (Classifiers[j] == 'LAmbDA'){
      res <- evaluate(truefile, predfile, Indices = c((Ind_SM2[i]+Train_Idx[1]):(Ind_SM2[i+1]+Train_Idx[1]-1)))
    } else if(Classifiers[j] == 'Garnett_CV'){
      res <- evaluate(truefile, predfile, Indices = c(Gar_SM2[i]:(Gar_SM2[i+1]-1)))
    } else {
      res <- evaluate(truefile, predfile, Indices = c(Ind_SM2[i]:(Ind_SM2[i+1]-1)))
    }
    
    MedF1_Big[j,Count_Big] <- res$MedF1
    Unl_Big[j,Count_Big] <- res$PercUnl
    Count_Big <- Count_Big + 1

  }
  
  
  for (i in c(1:(length(Ind_SW)-1))){
    truefile <- paste("results/PBMC_benchmarking/SW/",Classifiers[j],"/",Filename[j],"_True_Labels.csv",sep="")
    predfile <- paste("results/PBMC_benchmarking/SW/",Classifiers[j],"/",Filename[j],"_Pred_Labels.csv",sep="")
    
    if (Classifiers[j] == 'LAmbDA'){
      res <- evaluate(truefile, predfile, Indices = c((Ind_SW[i]+Train_Idx[7]):(Ind_SW[i+1]+Train_Idx[7]-1)))
    } else if(Classifiers[j] == 'Garnett_CV'){
      res <- evaluate(truefile, predfile, Indices = c(Gar_SW[i]:(Gar_SW[i+1]-1)))
    } else {
      res <- evaluate(truefile, predfile, Indices = c(Ind_SW[i]:(Ind_SW[i+1]-1)))
    }

      MedF1_Big[j,Count_Big] <- res$MedF1
      Unl_Big[j,Count_Big] <- res$PercUnl
      Count_Big <- Count_Big + 1

  }
  
  
}

#### Modify data for scoring ####
MedF1_Big <- gather_data2(MedF1_Big, Names, Train, Test, Exp_Big)
new_score = c(MedF1_Big$score,mean(c(0.90,0.99,0.99,0.96,0.84,0.90,0.99,
                                     0.99,0.96,0.97,0.96,0.98,0.95,
                                     0.96,0.40,0.63,0.45,0.95,0.92,0.99,
                                     1.0,0.92,0.99,0.97,0.67,0.87,0.99,
                                     1.0,0.99,0.99,0.96,0.96,0.91,0.98,
                                     1.0,1.0,0.93,1.0,1.0,1.0,1.0,
                                     1.0,0.99,1.0,0.96,0.95,1.0,0.98)))
names(new_score) <- c("SCALT",names(MedF1_Big$score))
#View(MedF1_Big$score)
#score <- order(new_score)
score <- seq(length(new_score),1)

Unl_Big <- gather_data2(Unl_Big, Names, Train, Test, Exp_Big)
Unl_Big <- Unl_Big$data

MedF1_Big$Unlab = Unl_Big$Value

#### Modify data for heatmap ####
MedF1_Big <- MedF1_Big$data
SCALT_results = head(MedF1_Big,48)
SCALT_results <- SCALT_results[,c(1,2,3)]

SCALT_results$Classifier <- rep("SCALT",48)
# SCALT_results$Value <- c(0.84,0.84,0.82,1.0,0.80,0.93,0.93,0.98,0.89,0.89,0.93,1.0,1.0,0.80,0.99,0.96,1.0,0.98,1.0,0.67,
#                          0.79,0.83,0.80,0.79,1.0,0.82,0.79,0.76,0.84,0.81,0.56,0.80,0.82,0.78,1.0,0.99,1.0,1.0,1.0,1.0,1.0,
#                          0.94,0.97,0.95,1.0,0.88,1.0,0.93)
# SCALT_results$Value <- c(0.87,0.86,0.79,0.88,0.86,0.92,0.90,0.89,0.90,0.84,0.89,0.99,0.87,0.38,1.0,0.67,0.92,
#                          0.94,0.88,0.11,0.82,0.70,0.58,0.70,0.83,0.74,0.86,0.97,0.97,0.95,0.55,0.96,0.96,0.95,1.0,0.93,
#                          0.68,1.0,1.0,1.0,0.72,0.92,0.92,0.91,0.85,0.93,0.94,0.96)
SCALT_results$Value <- c(0.90,0.99,0.99,0.96,0.84,0.90,0.99,
                         0.99,0.96,0.97,0.96,0.98,0.95,
                         0.96,0.40,0.63,0.45,0.95,0.92,0.99,
                         1.0,0.92,0.99,0.97,0.67,0.87,0.99,
                         1.0,0.99,0.99,0.96,0.96,0.91,0.98,
                         1.0,1.0,0.93,1.0,1.0,1.0,1.0,
                         1.0,0.99,1.0,0.96,0.95,1.0,0.98)

#MedF1_Big <- rbind(MedF1_Big,SCALT_results)
#New_Exp_Big <- c(Exp_Big,SCALT_results$Exp)
#New_Test <- c(Test,SCALT_results$Test)
#new_names <- c(Names,"SCALT")
#new_Expr <- c(Expr,"SCALT")

MedF1_Big <- rbind(SCALT_results,MedF1_Big)
New_Exp_Big <- c(SCALT_results$Exp,Exp_Big)
New_Test <- c(SCALT_results$Test,Test)
new_names <- c("SCALT",Names)
new_Expr <- c("SCALT",Expr)

factor_column = factor(MedF1_Big$Classifier,unique(MedF1_Big$Classifier))
MedF1_Big$block = factor_column

## heatmap -----
colors = mako(500)

p2 = ggplot(data = MedF1_Big, aes(x = Exp, y = block, fill = Value)) +
  geom_tile(colour="white",size=0.25) + geom_text(aes(label = round(Value,digits=2)), color = "Black", size = 2.5,fontface = "bold") +
  scale_fill_gradientn(colours=colors, limits=c(0, 1), guide = FALSE) +
  scale_y_discrete(expand=c(0,0), limits=new_names[score], labels = new_Expr[score]) + 
  scale_x_discrete(expand=c(0,0), breaks = New_Exp_Big, labels = New_Test, name = 'Test set') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.background=element_blank(),
        panel.border=element_blank()) +
  facet_grid( ~ Train, scales = "free_x", space = "free_x") + 
  geom_vline(xintercept = 6.51, colour = "red3", size = 1.1) +
  labs(fill = 'Median F1-score')
p2 = p2+scale_fill_bs5("green",reverse = T)
ggsave("significant_inter_pbmc_Abdelaal_benchmark.pdf",plot =p2,dpi = 1000,width = 14, height = 6)

#### Proportion of unlabeled ###
unlabScalt <- head(Unl_Big,48)
unlabScalt <- unlabScalt[,c(1,2,3)]
unlabScalt$Classifier <- rep("SCALT",48)
# unlabScalt$Value <- c(0.0,0.10,0.02,0.68,0.87,0.84,0.0,0.01,0.49,0.06,0.75,0.91,0.88,0.15,1.0,0.99,1.0,0.99,
#                       0.97,0.01,0.04,0.69,0.42,0.05,0.94,0.86,0.28,0.36,0.74,0.54,0.14,0.88,0.88,0.83,1.0,1.0,
#                       0.21,1.0,1.0,1.0,0.0,0.01,0.59,0.20,0.03,0.59,0.89,0.85)
unlabScalt$Value <- c(0.0,0.12,0.0,0.55,0.85,0.74,0.25,0.0,0.37,0.0,0.61,0.91,0.82,0.12,0.97,0.85,0.97,0.99,0.99,0.0,
                      0.0,0.41,0.10,0.0,0.89,0.71,0.16,0.21,0.52,0.35,0.0,0.66,0.77,0.75,0.98,0.93,0.17,0.99,0.99,0.99,0.0,
                      0.0,0.27,0.0,0.0,0.29,0.80,0.70)

Unl_Big <- rbind(unlabScalt,Unl_Big)
factor_column = factor(Unl_Big$Classifier,unique(Unl_Big$Classifier))
Unl_Big$block = factor_column

colors = cividis(500) 
p4 = ggplot(data = Unl_Big, aes(x = Exp, y = block, fill = Value)) +
  geom_tile(colour="white",size=0.25) + geom_text(aes(label = round(Value,digits=2)), color = "Black", size = 2.0,fontface = "bold") +
  scale_fill_gradientn(colours=colors, limits=c(0, 1), guide = FALSE) +
  scale_y_discrete(expand=c(0,0), limits=new_names[score], labels = new_Expr[score]) +
  scale_x_discrete(expand=c(0,0), breaks = New_Exp_Big, labels = New_Test, name = 'Test set') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.background=element_blank(),
        panel.border=element_blank()) +
  facet_grid( ~ Train, scales = "free_x", space = "free_x") +
  geom_vline(xintercept = 6.5, color = 'red3', size = 1.1) +
  labs(fill = 'Unlabeled (%)')
p4 = p4 + scale_fill_bs5("indigo")
ggsave("significant_inter_pbmc_unlabelledAbdelaal_benchmark.pdf",plot =p4,dpi = 1000,width = 14, height = 6)

## boxplots -----

p <- ggplot(MedF1_Big, aes(x = Exp, y = Value)) + 
  geom_boxplot() + 
  xlab("Dataset") + 
  ylab("Median F1-score") + 
  scale_x_discrete(expand=c(0,0), breaks = New_Exp_Big, labels = New_Test) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = 'black'),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  facet_grid( ~ Train, scales = "free_x", space = "free_x") +
  geom_vline(xintercept = 6.5, color = "red3", size = 1.1)

ggsave("significant_x_boxplots.pdf",plot =p,dpi = 1000,width = 13, height = 2.3)

p <- ggplot(MedF1_Big, aes(x = Classifier, y = Value)) + 
  geom_boxplot() + 
  stat_summary(fun = 'mean', colour = "red", size = 2, geom = 'point') +
  xlab("Classifier") + 
  ylab("Median F1-score") + 
  scale_x_discrete(expand=c(0,0), limits = new_names[score]) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), labels = c('0','0.25','0.5','0.75','1'))+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.line = element_line(colour = 'black'),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  
  coord_flip()
ggsave("significant_y_boxplots.pdf",plot =p,dpi = 1000,width = 2, height = 6)

# #### La mediana è biased ####
# Train_labs_biased <- c('10XV2', '10XV2', '10XV2', '10XV2', '10XV2', '10XV2', '10XV2',
#            '10XV3', '10XV3', '10XV3', '10XV3', '10XV3', '10XV3',
#            'CL', 'CL', 'CL', 'CL', 'CL', 'CL', 'CL',
#            'DR', 'DR', 'DR', 'DR', 'DR', 'DR', 'DR',
#            'iD', 'iD', 'iD', 'iD', 'iD', 'iD', 'iD',
#            'SM2', 'SM2', 'SM2', 'SM2', 'SM2', 'SM2', 'SM2',
#            'SW', 'SW', 'SW', 'SW', 'SW', 'SW', 'SW')
# 
# Test_labs_biased <- c('SM2', '10XV3', 'CL', 'DR', 'iD', 'SW', '10XV2',
#           'SM2', '10XV2', 'CL', 'DR', 'iD', 'SW',
#           'SM2', '10XV2', '10XV3', 'DR', 'iD', 'SW', 'CL',
#           'SM2', '10XV2', '10XV3', 'CL', 'iD', 'SW', 'DR',
#           'SM2', '10XV2', '10XV3', 'CL', 'DR', 'SW', 'iD',
#           '10XV2', '10XV3', 'CL', 'DR', 'iD', 'SW', 'SM2',
#           'SM2', '10XV2', '10XV3', 'CL', 'DR', 'iD', 'SW')
# 
# prop_zeros = c(0/6,0/8,0/7,1/8,2/7,0/7,0/9,0/6,0/8,0/6,1/8,2/6,0/6,0/6,4/7,2/7,5/7,4/6,2/4,0/5,1/6,1/9,1/8,0/7,3/7,1/7,0/6,
#          1/6,1/7,0/6,0/6,2/6,1/6,0/7,4/6,4/6,1/6,4/6,6/7,2/5,1/5,0/5,0/7,0/6,0/5,0/7,1/6,0/4)
# median_SCALT = c(0.82,0.85,0.79,0.80,0.49,0.73,0.77,0.89,0.87,0.81,0.83,0.82,0.89,0.80,0.09,0.13,0.0,0.0,0.6,0.52,
#                  0.81,0.67,0.69,0.67,0.22,0.72,0.83,0.92,0.87,0.87,0.78,0.94,0.89,0.87,0.0,0.0,0.67,0.0,0.0,1.0,0.83,
#                  0.93,0.97,0.90,0.84,0.85,0.68,0.92)
# unlabelled_SCALT =  c(0.0,0.12,0.0,0.55,0.85,0.74,0.25,0.0,0.37,0.0,0.61,0.91,0.82,0.12,0.97,0.85,0.97,0.99,0.99,0.0,
#                                          0.0,0.41,0.10,0.0,0.89,0.71,0.16,0.21,0.52,0.35,0.0,0.66,0.77,0.75,0.98,0.93,0.17,0.99,0.99,0.99,0.0,
#                                          0.0,0.27,0.0,0.0,0.29,0.80,0.70)
# Exp_Big_biased <- c('10XV2 - SM2', '10XV2 - 10XV3', '10XV2 - CL', '10XV2 - DR', '10XV2 - iD', '10XV2 - SW', '10XV2 - 10XV2',
#              '10XV3 - SM2', '10XV3 - 10XV2', '10XV3 - CL', '10XV3 - DR', '10XV3 - iD', '10XV3 - SW',
#              'CL - SM2', 'CL - 10XV2', 'CL - 10XV3', 'CL - DR', 'CL - iD', 'CL - SW', 'CL - CL',
#              'DR - SM2', 'DR - 10XV2', 'DR - 10XV3', 'DR - CL', 'DR - iD', 'DR - SW', 'DR - DR',
#              'iD - SM2', 'iD - 10XV2', 'iD - 10XV3', 'iD - CL', 'iD - DR', 'iD - SW', 'iD - iD',
#              'SM2 - 10XV2', 'SM2 - 10XV3', 'SM2 - CL', 'SM2 - DR', 'SM2 - iD', 'SM2 - SW', 'SM2 - SM2',
#              'SW - SM2', 'SW - 10XV2', 'SW - 10XV3', 'SW - CL', 'SW - DR', 'SW - iD', 'SW - SW')
# 
# biased_median = data.frame(Value = c(median_SCALT,unlabelled_SCALT,prop_zeros),
#                            Class=c(rep("Median F1",48),rep("% Unlabeled",48),rep("% Zeros",48)),
#                            Exp = rep(Exp_Big_biased,3),Train = rep(Train_labs_biased,3),Test=rep(Test_labs_biased,3))
# 
# 
# colors = plasma(500)
# p5 = ggplot(data = biased_median, aes(x = Exp, y = Class, fill = Value))+geom_tile(colour="white",size=0.25)+
#   scale_fill_gradientn(colours=colors, limits=c(0, 1),guide = FALSE)+ 
#   geom_text(aes(label = round(Value,digits=2)), color = "Black", size = 2.8,fontface = "bold")+
#   scale_x_discrete(labels= New_Test)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         strip.text.x = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         plot.background=element_blank(),
#         panel.border=element_blank(),aspect.ratio = 1)+
#   facet_grid( ~ Train,scales = "free_x",)+#,  space = "free_x"
#   geom_vline(xintercept = 6.5, color = 'red3', size = 1.1) +
#   labs(fill = 'Unlabeled (%)')
# 
# ggsave("median_problem.pdf",plot =p5,dpi = 1000,width = 14, height = 8)
