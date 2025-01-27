#!/usr/bin/env Rscript
library(RColorBrewer)
library(stringr) 
library(ggplot2)
library(viridis)

barplot_function <- function(ass,multi,unasss){
  A = data.frame(table(ass$Expression))
  colnames(A)=c("Class","Value")
  M = data.frame(table(multi$Expression))
  colnames(M)=c("Class","Value")
  U = data.frame(table(unasss$Expression))
  colnames(U)=c("Class","Value")
  Class = c(rep("Assigned",2),rep("Multiassigned",2),rep("Unassigned",2))
  Condition = rep(c("Genes expressed < 500","Genes expressed >= 500"),3)
  Values = c(A$Value,M$Value,U$Value)
  data = data.frame(Class=Class,Condition=Condition,Values=Values)
  bar_plot = ggplot(data, aes(fill=Condition, y=Values, x=Class)) + ggtitle("")+
    geom_bar(position="dodge", stat="identity")+scale_fill_brewer(palette = "Dark2")+theme_classic()+
    theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.title = element_blank())
  return(bar_plot)
}

df = read.table("../last_merge/config_index.tsv",header = F,sep="\t")
colnames(df)=c("Original","SCALT","Cell-ID","Sample","Expression")

#### Fibroblast barplot ####
fibro = df[grepl("FIBROBLAST",df[,"Original"]),]
multiassigned = fibro[grepl("multiassigned",fibro[,"SCALT"]),]
unssigned = fibro[grepl("unassigned",fibro[,"SCALT"]),]
assigned = fibro[!grepl("multiassigned|unassigned",fibro[,"SCALT"]),]
fibro_barplot = barplot_function(assigned,multiassigned,unssigned)
png("fibro_barplot.png",width= 10,height= 6,units= "in",res=1200)
fibro_barplot
dev.off()

#### Neuron barplot ####
neuro = df[grepl("NEURON|GLIA|ASTROCYTE",df[,"Original"]),]
multiassigned = neuro[grepl("multiassigned",neuro[,"SCALT"]),]
unssigned = neuro[grepl("unassigned",neuro[,"SCALT"]),]
assigned = neuro[!grepl("multiassigned|unassigned",neuro[,"SCALT"]),]
neuro_barplot = barplot_function(assigned,multiassigned,unssigned)
png("neuro_barplot.png",width= 10,height= 6,units= "in",res=1200)
neuro_barplot
dev.off()

#### Immune barplots ####
tcells_notmegred = c("DOUBLE.NEGATIVE.T.CELL","T.CELL")
bcells_notmerged = c("B-CELLS.HPA.CELL","PREB.CELL","PROB.CELL")
myeloid_notmerged = c("CDC1.CELL","COMMON.LYMPHOID.PROGENITOR.CELL","COMMON.MYELOID.PROGENITOR.CELL","CYCLING.CDC2.CELL"
                      ,"DENDRITIC.CELLS.HPA.CELL","LANGERHANS.CELLS.HPA.CELL","MACROPHAGES.HPA.CELL",
                      "MAST.CELL","MEGAKARYOCYTE.CELL","MEGAKARYOCYTE.PROGENITOR.CELL","MHC.CLASS.II.HIGH.MONOCYTE.DC.MACROPHAGE.CELL",
                      "PLACENTA.CAMP+.ANTIMICROBIAL.MONOCYTE.CELL","PLACENTA.DEFENSIN+.ANTIMICROBIAL.MONOCYTE.CELL","PRE-PDC..LYMPHOID.ORIGIN..CELL",
                      "PROLIFERATION.MAST.CELL","PDC.CELL","ITGAX+.MACROPHAGE.CELL","CYCLING.MYELOID.CELL", "PRE-PDC.MYELOID.CELL",
                      "PROLIFERATION.MYELOID.CELL","PLACENTA.CAMP+.ANTIMICROBIAL.MONOCYTE.CELL",
                      "PLACENTA.DEFENSIN+.ANTIMICROBIAL.MONOCYTE.CELL","ITGAX+.MACROPHAGE.CELL")
neutrophils_notmerged = c("METAMYELOCYTE.BAND.NEUTROPHIL.CELL")
hofbauer_notmerged = c("HOFBAUER.CELL")
blymph_pro_notmerged = c("LYMPHOID.PRE-PDC.CELL")

merged = c("MULTIFUNCTIONAL.T.CELL.2","MULTIFUNCTIONAL.B.CELL","MYELOID.CELL.2","NEUTROPHIL.GRANULOGYTE.CELL",
           "ALVEOLAR.MACROPHAGE.CELL","TISSUE.RESIDENT.MACROPHAGE.CELL","PROLIFERATION.HOFBAUER.CELL",
           "B.LYMPHOCYTE.PROGENITOR.CYCLING.CELL","CYCLING.DN.DP.T.CELL","DARK.LIGHT.ZONE.GC.B.CELL","BREAST.MILK.MACROPHAGE.CELL",
           "DOUBLE.POSITIVE.CELL","PROLIFERATION.DIFFERENTING.T.CELL","CYCLING.MYELOID.CELL","PRE-PDC.MYELOID.CELL","PROLIFERATION.MYELOID.CELL")

ordered_cells_to_match = c(merged,tcells_notmegred,bcells_notmerged,myeloid_notmerged,neutrophils_notmerged,
                           hofbauer_notmerged,blymph_pro_notmerged)

nnnn = rep(NA,length(ordered_cells_to_match))
for (y in 1:length(nnnn)){
  nnnn[y]=paste(paste("\\<",ordered_cells_to_match[y],sep=""),"\\>",sep="")
}

togrep = paste(nnnn,collapse = "|")

immune = df[grepl(togrep,df[,"Original"]),]
multiassigned = immune[grepl("multiassigned",immune[,"SCALT"]),]
unssigned = immune[grepl("unassigned",immune[,"SCALT"]),]
assigned = immune[!grepl("multiassigned|unassigned",immune[,"SCALT"]),]
immune_barplot = barplot_function(assigned,multiassigned,unssigned)
png("immune_barplot.png",width= 10,height= 6,units= "in",res=1200)
immune_barplot
dev.off()



