#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("pheatmap")
library(RColorBrewer)
library(stringr) 
library("ComplexHeatmap")
library("ggsci")

confusion_matrix <-function(TrueLabelsPath,PredLabelsPath){
  Indices = NULL
  true_lab <- unlist(read.csv(TrueLabelsPath))
  pred_lab <- unlist(read.csv(PredLabelsPath))
  
  if (! is.null(Indices)){
    true_lab <- true_lab[Indices]
    pred_lab <- pred_lab[Indices]
  }
  
  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))
  
  unique_all <- unique(c(unique_true,unique_pred))
  conf <- as.data.frame.matrix(table(true_lab,pred_lab))
  scaled_conf = conf/apply(conf,1,sum)
  
  new_row_names = str_to_title(tolower(rownames(scaled_conf)))
  new_col_names = str_to_title(tolower(colnames(scaled_conf)))
  return(scaled_conf)
}

barplot_matrix <- function(TrueLabelsPath,PredLabelsPath){
  Indices = NULL
  true_lab <- unlist(read.csv(TrueLabelsPath))
  pred_lab <- unlist(read.csv(PredLabelsPath))
  
  if (! is.null(Indices)){
    true_lab <- true_lab[Indices]
    pred_lab <- pred_lab[Indices]
  }
  
  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))
  
  unique_all <- unique(c(unique_true,unique_pred))
  conf <- as.data.frame.matrix(table(true_lab,pred_lab))
  sums = data.frame("Rel. Aboundande"=apply(conf,1,sum))
  return(sums)
}

mypal <- pal_npg("nrc", alpha = 1.0)(9)
#scales::show_col(mypal)
# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF"

#TrueLabelsPath = args[1]
#PredLabelsPath = args[2]
TrueLabelsPath = "base_case/try_01/original.csv"
PredLabelsPath = "base_case/try_01/annotated.csv"
scaled_conf=confusion_matrix(TrueLabelsPath,PredLabelsPath)

#### Fibroblast with initial labels and without multiassigned and unassigned####
fibroblast_conf = scaled_conf[grepl("FIBROBLAST",rownames(scaled_conf)),grepl("FIBROBLAST",colnames(scaled_conf))]
fibro_new_row_names = str_to_title(tolower(rownames(fibroblast_conf)))
fibro_new_col_names = str_to_title(tolower(colnames(fibroblast_conf)))
rownames(fibroblast_conf)=fibro_new_row_names
colnames(fibroblast_conf)=fibro_new_col_names

o = c("Postn.fibroblast.cell","Nr4a3.Fibroblast.cell","Scn7a.fibroblast.cell","Pcolce2.Fibroblast.cell",
      "Metallothionein+.Fibroblast.cell","Tnc.fibroblast.cell","Igfbp6+Apod+.Fibroblast.cell","Clu+.Fibroblast.cell",
      "Cfd.fibroblast.cell","Fibroblasts.hpa.cell","Ptgds+Cxcl14+.Fibroblast.cell","Thy1+.Fibroblast.cell","Wisp2.Fibroblast.cell",
      "Cd55+Sema3c+.Fibroblast.cell","Igfbp6+Sfrp4+.Fibroblast.cell","Cycling.fibroblast.cell","C1qtnf3+Shisa3+Col14a1+.Fibroblast.cell",
      "C7+.Fibroblast.cell","Cxcl14.Fibroblast.cell","Proliferation.cxcl14.Fibroblast.cell")

new_names = c(o,setdiff(fibro_new_row_names,o))
ordered_fibroblast_conf = fibroblast_conf[new_names,new_names]

adj_name_fibro = c("POSTN fibroblast","NR4A3 Fibroblast","SCN7A fibroblast","PCOLCE2 Fibroblast",
                   "Metallothionein+ Fibroblast","TNC fibroblast","IGFBP6+ APOD+ Fibroblast",
                   "CLU+ Fibroblast","CFD fibroblast","Fibroblasts (HPA)","PTGDS+ CXCL14+ Fibroblast",
                   "THY1+ Fibroblast","WISP2 Fibroblast","CD55+ SEMA3C+ Fibroblast",
                   "IGFBP6+ SFRP4+ Fibroblast","Cycling fibroblast","C1QTNF3+ SHISA3+ COL14A1+ Fibroblast",
                   "C7+ Fibroblast","CXCL14 Fibroblast","Proliferation CXCL14 Fibroblast",
                   "ADAMDEC1+ CXCL14+ Fibroblast","Adipose ITLN1 Fibroblast","ANGPTL7+ NRP2+ Fibroblast",
                   "APCDD1 Fibroblast","APOD fibroblast","APOE fibroblast","C2orf40 Fibroblast",
                   "CCL19.21 Fibroblast","CCL2+ SFRP2+ Fibroblast","CDCA7 Proliferation fibroblast",
                   "COCH fibroblast","COL11A1 Fibroblast","CXCL+ Fibroblast","CXCL1.2.3 Fibroblast",
                   "CXCL13 Fibroblast","CXCL14+ POSTN+ Fibroblast","F13A1 Fibroblast",
                   "FBLN2+ APOD+ Fibroblast","FGFBP2 Fibroblast","Fibroblast cell",
                   "Fibroblast granulosa doublet like cell","FRZB fibroblast cell",
                   "G0S2+ PPP1R14A+ Fibroblast","GREB1 Fibroblast","IGFBP1 Fibroblast",
                   "Metallothionein+ Collagen- Fibroblast","MKI67 Proliferation fibroblast",
                   "MMP1 Fibroblast cell","Myofibroblast","PDGFRA+ BMP4+ Fibroblast",
                   "PLA2G2A+ IGFBP4+ APOC1+ Fibroblast","PLAC9 Fibroblast","RNASE1 Fibroblast",
                   "SPARCL1 Fibroblast","TCIM+ NRG1+ Fibroblast","THY1+ Fibroblast Like pericyte")  
row.names(ordered_fibroblast_conf)=adj_name_fibro
colnames(ordered_fibroblast_conf)=adj_name_fibro
## BARS ## #8c6bb1
bars = barplot_matrix(TrueLabelsPath,PredLabelsPath)
rownames(bars)=str_to_title(tolower(rownames(bars)))
filt_bars = log10(data.frame("Rel. Aboundande"=bars[new_names,]))
ha_list = rowAnnotation(" "= anno_barplot(filt_bars, gp = gpar(fill = "#8491B4FF"), #"#084081"
                                          bar_width = 1, width = unit(2, "cm")))

ordered_fibroblast_conf[ordered_fibroblast_conf==0] <- NA

## HEATMAP ## "#ccebc5","#7bccc4","#0868ac"
png("fibro_initial_labels_plot.png",width= 10,height= 6,units= "in",res=1200)
PLOT1 = Heatmap(ordered_fibroblast_conf,cluster_columns = F,cluster_rows = F,col=colorRampPalette(c("#E8EAF6","#5C6BC0","#303F9F"))(100),
        column_names_rot = 45,na_col = "#f7fcf0",row_names_gp = gpar(fontsize = 5,fontface = "bold"),column_names_gp = gpar(fontsize = 5,fontface = "bold"),
        heatmap_legend_param = list(title = "",legend_direction = "horizontal",legend_width = unit(5, "cm")),row_names_side = c("left"),
        right_annotation = ha_list)
draw(PLOT1, heatmap_legend_side = "bottom")
dev.off()

rm(list=setdiff(ls(), c("confusion_matrix","barplot_matrix")))

#### Fibroblasts after merging and with Multiassigned, assigned and unassigend categories ####

TrueLabelsPath = "last_merge/try_02/original.csv"
PredLabelsPath = "last_merge/try_02/annotated.csv"
scaled_conf=confusion_matrix(TrueLabelsPath,PredLabelsPath)

fibroblast_conf = scaled_conf[grepl("FIBROBLAST",rownames(scaled_conf)),grepl("FIBROBLAST",colnames(scaled_conf))]
fibro_new_row_names = str_to_title(tolower(rownames(fibroblast_conf)))
fibro_new_col_names = str_to_title(tolower(colnames(fibroblast_conf)))
rownames(fibroblast_conf)=fibro_new_row_names  
colnames(fibroblast_conf)=fibro_new_col_names

adj_names_fibro_pruned = c("ADAMDEC1+ CXCL14+ Fibroblast","ADIPOSE ITLN1 Fibroblast",
                           "ANGPTL7+ NRP2+ Fibroblast","APCDD1 Fibroblast","APOD fibroblast",
                           "APOE fibroblast","C2orf40 Fibroblast","CCL19.21 Fibroblast",
                           "CCL2+ SFRP2+ Fibroblast","COCH fibroblast","COL11A1 Fibroblast",
                           "CXCL+ Fibroblast","CXCL1.2.3 Fibroblast","CXCL13 Fibroblast",
                           "CXCL14 Fibroblast","CXCL14+ POSTN+ Fibroblast","F13A1 Fibroblast",
                           "FBLN2+ APOD+ Fibroblast","FGFBP2 Fibroblast","Fibroblast cell",
                           "FRZB fibroblast","G0S2+ PPP1R14A+ Fibroblast","General Fibroblast",
                           "GREB1 Fibroblast","Heart Fibroblast","IGFBP1 Fibroblast",
                           "Intestinal Fibroblast","Metallothionein+ Collagen- Fibroblast","MMP1 Fibroblast",
                           "Myofibroblast","PDGFRA+ BMP4+ Fibroblast","PLA2G2A+ IGFBP4+ APOC1+ Fibroblast",
                           "SPARCL1 Fibroblast","TCIM+ NRG1+ Fibroblast","THYM1+ Fibroblast Like pericyte")  

colnames(fibroblast_conf)=adj_names_fibro_pruned
row.names(fibroblast_conf)=adj_names_fibro_pruned

toMatch <- c("FIBROBLAST","unassigned","multiassigned")
m4barplot = scaled_conf[grepl(paste(toMatch,collapse="|"),rownames(scaled_conf)),grepl(paste(toMatch,collapse="|"),colnames(scaled_conf))]
bars = data.frame("Assigned"=1-apply(m4barplot[,c("multiassigned","unassigned")],1,sum),
                  "Multiassigned"= m4barplot["multiassigned"],"Unassigned"=m4barplot["unassigned"])

ha_list = rowAnnotation(" "= anno_barplot(bars, gp = gpar(fill = c("#fdbf6f","#cab2d6","#b2df8a")), #"#f7fcf0","#7bccc4","#0868ac"
                                          bar_width = 1, width = unit(2, "cm")))

png("cleaned_lists_merged_fibro_plot.png",width= 10,height= 6,units= "in",res=1200)
PLOT2 <- Heatmap(fibroblast_conf,cluster_columns = T,cluster_rows = T,show_row_dend = FALSE,
                show_column_dend = FALSE,row_names_gp = gpar(fontsize = 6,fontface = "bold"),column_names_gp = gpar(fontsize = 6,fontface = "bold"),column_names_side = c("bottom"),
                row_names_side = c("left"),right_annotation = ha_list,column_names_rot = 45,na_col = "#f7fcf0",
                heatmap_legend_param = list(title = ""),col=colorRampPalette(c("#E8EAF6","#5C6BC0","#303F9F"))(100))#col = colorRamp2(c(0, 0.5, 1),c("#f7f7f7","#c2a5cf","#762a83"))
lgd = Legend(labels = c("Assigned","Multiassigned","Unassigned"), title = "Class", legend_gp = gpar(fill = c("#fdbf6f","#cab2d6","#b2df8a")))
draw(PLOT2,annotation_legend_list = lgd) 
dev.off()

rm(list=setdiff(ls(), c("confusion_matrix","barplot_matrix")))

#### NEURONS - initial assignement without multiassigned and unassigned ####
TrueLabelsPath = "base_case/try_01/original.csv"
PredLabelsPath = "base_case/try_01/annotated.csv"
scaled_conf=confusion_matrix(TrueLabelsPath,PredLabelsPath)

toMatch <- c("NEURON","GLIA","ASTROCYTE")
neuron_conf = scaled_conf[grepl(paste(toMatch,collapse="|"),rownames(scaled_conf)),grepl(paste(toMatch,collapse="|"),colnames(scaled_conf))]
neuro_new_row_names = str_to_title(tolower(rownames(neuron_conf)))
neuro_new_col_names = str_to_title(tolower(colnames(neuron_conf)))
rownames(neuron_conf)=neuro_new_row_names
colnames(neuron_conf)=neuro_new_col_names
rownames(neuron_conf)  #1:9;10:13;14:16;17:19,20:21,22:23
o = c("L5.6.Excitatory.neuron.cell","Plch1.L4.5.Excitatory.neuron.cell","L4.Excitatory.neuron.cell",
      "L2.3.Excitatory.neuron.cell","L5b.excitatory.neuron.cell","L5.Excitatory.neuron.cell",
      "Tshz2.L4.5.Excitatory.neuron.cell","L6.Excitatory.neuron.cell","Pyramidal.neuron.cell",
      "Glia.progenitor.cell","S.phase.glial.cell","G2.M.phase.glial.cell","Apoe+Bcan+.Glial.cell",
      "Cxcl14.Inhibitory.neuron.cell","Sv2c.inhibitory.neuron.cell","Vip.inhibitory.neuron.cell",
      "Glial.cell","Mbp+Mal+Aatk+.Glial.cell","Cryab+Plat+.Glial.cell","Fibrous.astrocyte.cell",
      "Protoplasmic.astrocyte.cell","Bnc2+.Neuron.cell","Etv1+.Neuron.cell")

new_names = c(o,setdiff(neuro_new_row_names,o))
ordered_neuro_conf = neuron_conf[new_names,new_names]
row.names(ordered_neuro_conf)

adj_names = c("L5.6 Excitatory neuron","PLCH1 L4.5 Excitatory neuron","L4 Excitatory neuron",
              "L2.3 Excitatory neuron","L5b Excitatory neuron","L5 Excitatory neuron",
              "TSHZ2 L4.5 Excitatory neuron","L6 Excitatory neuron","Pyramidal neuron","Glia progenitor",
              "S phase glial","G2M phase glial","APOE+BCAN+ Glial","CXCL14 Inhibitory neuron",
              "SV2C Inhibitory neuron","VIP Inhibitory neuron","Glial cell","MBP+MAL+AATK+ Glial cell",
              "CRYAB+PLAT+ Glial cell","Fibrous astrocyte","Protoplasmic astrocyte","BNC2+ Neuron",
              "ETV1+ Neuron","Astrocyte","Astrocytes (HPA)","Excitatory neurons (HPA)",
              "Inhibitory neurons (HPA)","Microglia cell","Microglial cells (HPA)","Muller glia",
              "Muller glia.cells (HPA)","Neuron","PVALB Inhibitory neuron","SST Inhibitory neuron") 
row.names(ordered_neuro_conf)=adj_names
colnames(ordered_neuro_conf)=adj_names

## BARS ## #023858
bars = barplot_matrix(TrueLabelsPath,PredLabelsPath)
rownames(bars)=str_to_title(tolower(rownames(bars)))
filt_bars = log10(data.frame("Rel. Aboundande"=bars[new_names,]))
ha_list = rowAnnotation(" "= anno_barplot(filt_bars, gp = gpar(fill = "#4DBBD5FF"),
                                          bar_width = 1, width = unit(2, "cm")))

## HEATMAP ## "#ece7f2","#3690c0","#045a8d"
ordered_neuro_conf[ordered_neuro_conf==0] <- NA
png("neuro_initial_labels_plot.png",width= 10,height= 6,units= "in",res=1200)
PLOT1 = Heatmap(ordered_neuro_conf,cluster_columns = F,cluster_rows = F,col=colorRampPalette(c("#E3F2FD","#42A5F5","#1976D2"))(100),
        column_names_rot = 45,na_col = "#fff7fb",row_names_gp = gpar(fontsize = 7,fontface = "bold"),column_names_gp = gpar(fontsize = 7,fontface = "bold"),
        heatmap_legend_param = list(title = "",legend_direction = "horizontal",legend_width = unit(5, "cm")),row_names_side = c("left"),
        right_annotation = ha_list)
draw(PLOT1, heatmap_legend_side = "bottom")
dev.off()

rm(list=setdiff(ls(), c("confusion_matrix","barplot_matrix")))

#### Neurons after merging and with Multiassigned, assigned and unassigend categories ####
TrueLabelsPath = "last_merge/try_02/original.csv"
PredLabelsPath = "last_merge/try_02/annotated.csv"
scaled_conf=confusion_matrix(TrueLabelsPath,PredLabelsPath)

neuroToMatch = toMatch <- c("NEURON","GLIA","ASTROCYTE")
neuro_conf = scaled_conf[grepl(paste(neuroToMatch,collapse="|"),rownames(scaled_conf)),grepl(paste(neuroToMatch,collapse="|"),colnames(scaled_conf))]
neuro_new_row_names = str_to_title(tolower(rownames(neuro_conf)))
neuro_new_col_names = str_to_title(tolower(colnames(neuro_conf)))
rownames(neuro_conf)=neuro_new_row_names  
colnames(neuro_conf)=neuro_new_col_names

toMatch <- c("NEURON","GLIA","ASTROCYTE","unassigned","multiassigned")
m4barplot = scaled_conf[grepl(paste(toMatch,collapse="|"),rownames(scaled_conf)),grepl(paste(toMatch,collapse="|"),colnames(scaled_conf))]
bars = data.frame("Assigned"=1-apply(m4barplot[,c("multiassigned","unassigned")],1,sum),
                  "Multiassigned"= m4barplot["multiassigned"],"Unassigned"=m4barplot["unassigned"])

ha_list = rowAnnotation(" "= anno_barplot(bars, gp = gpar(fill = c("#fdbf6f","#cab2d6","#b2df8a")), 
                                          bar_width = 1, width = unit(2, "cm")))
adj_names_pruned = c("Astrocyte","Astrocytes (HPA)","Excitatory neuron","Excitatory neurons (HPA)",
                     "Fibrous protoplasmatic astrocyte","General glial","Glial progenitor cycling",
                     "Inhibitory neuron","Inhibitory neurons (HPA)","Intestinal neuron","Microglia cell",
                     "Microglial cells (HPA)","Muller glia","Neuron","PVALB Inhibitory neuron",
                     "SST Inhibitory neuron")   
row.names(neuro_conf)=adj_names_pruned
colnames(neuro_conf)=adj_names_pruned
png("cleaned_lists_merged_neuro_plot.png",width= 10,height= 6,units= "in",res=1200)
PLOT2 <- Heatmap(neuro_conf,cluster_columns = F,cluster_rows = F,show_row_dend = FALSE,
                 show_column_dend = FALSE,row_names_gp = gpar(fontsize = 8,fontface = "bold"),column_names_gp = gpar(fontsize = 8,fontface = "bold"),column_names_side = c("bottom"),
                 row_names_side = c("left"),right_annotation = ha_list,column_names_rot = 45,na_col = "#f7fcf0",
                 heatmap_legend_param = list(title = ""),col=colorRampPalette(c("#E3F2FD","#42A5F5","#1976D2"))(100))#col = colorRamp2(c(0, 0.5, 1),c("#f7f7f7","#c2a5cf","#762a83"))
lgd = Legend(labels = c("Assigned","Multiassigned","Unassigned"), title = "Class", legend_gp = gpar(fill = c("#fdbf6f","#cab2d6","#b2df8a")))
draw(PLOT2,annotation_legend_list = lgd) 
dev.off()

rm(list=setdiff(ls(), c("confusion_matrix","barplot_matrix")))

#### IMMUNE SYSTEM ####
TrueLabelsPath = "base_case/try_01/original.csv"
PredLabelsPath = "base_case/try_01/annotated.csv"
scaled_conf=confusion_matrix(TrueLabelsPath,PredLabelsPath)

## T cells ##
tcells_merged = c("PROLIFERATION.IMMUNE.CELL","GAMMA.DELTA.T.CELL","CD16.NK.CELL","NK.CELL","CYCLING.T.CELL","GZMB.CD8.T.CELL",
           "CYTOTOXIC.CD8.T.CELL","CYCLING.T.NK.CELL","CD56.NK.CELL","GZMK+IL7R-.CD8.T.CELL","GZMK.NK.CELL","TISSUE-RESIDENT.NK.CELL",
           "NK-CELLS.HPA.CELL","T.NK.CELL","MEMORY.T.CELL","GZMK+IL7R+.CD8.T.CELL","KLRB1.CD8.T.CELL","CD8.T.CELL",
           "MEMORY.CD8.T.CELL","GZMK.CD8.T.CELL","IFN-ACTIVATED.T.CELL","PROLIFERATION.T.NK.CELL","CXCL13.EXHAUSTED.CD8.T.CELL",
           "CD4.T.CELL","TREG.CELL","KLRB1.CYTOTOXIC.CD4.T.CELL","GZMK.CYTOTOXIC.CD4.T.CELL","EFFECTOR.MEMORY.CD4.T.CELL",
           "MEMORY.CD4.T.CELL","NAIVE.CD4.T.CELL","NAIVE.T.CELL","INF.RESPONSED.T.CELL","DIFFERENTIATING.TREG.CELL",
           "NAIVE.CD8.T.CELL","AGONIST.T.CELL","CD8AA+.T.CELL","GDT.CELL","MAIT.CELL","T-CELLS.HPA.CELL")

tcells_notmegred = c("DOUBLE.NEGATIVE.T.CELL","T.CELL")

## B cells ##
bcells_merged = c("PRE-GC.B.CELL","GC.COMMITTED.CELL","INF.RESPONSED.NAIVE.B.CELL","NON-CLASS.SWITCH.MEMORY.B.CELL",
                  "DZ.LZ.B.CELL","CLASS.SWITCH.MEMORY.B.CELL","MEMORY.B.PRECURSOR.CELL","NAIVE.B.CELL","LZ.B.CELL",
                  "FCRL4.MEMORY.B.CELL","S.PHASE.DZ.B.CELL","G2M.DZ.B.CELL","HISTONE.HIGH.S.PHASE.DZ.B.CELL",
                  "NONPROLIFERATIVE.GC.B.CELL","PLASMA.CELL.PRECURSOR.CELL","PIF1+.G2M.DZ.B.CELL","B.CELL",
                  "MEMORY.B.CELL")
bcells_notmerged = c("B-CELLS.HPA.CELL","PREB.CELL","PROB.CELL")

## Myeloid ##
myeloid_merged = c("MONOCYTE.CELL","CD14.MONOCYTE.CELL","MHCII.HIGH.CD14.MONOCYTE.CELL","MHCII.LOW.CD14.MONOCYTE.CELL",
                   "CD16.MONOCYTE.CELL","CD14.CD16.MONOCYTE.CELL","CDC2.CELL","CD14+MHCIIHIGH.MONOCYTE.CELL","CD14+MHCIILOW.MONOCYTE.CELL",
                   "DENDRITIC.CELL","MACROPHAGE.CELL","S100A+.PRENEUTROPHIL..CYCLING..CELL","MYELOCYTE.CELL",
                   "CYCLING.S100A+.PRENEUTROPHIL.CELL","S100A+.PRENEUTROPHIL.CELL","PROMYELOCYTE.CELL","MONOCYTES.HPA.CELL",
                   "CDP.CELL","GMP.CELL","GRANULOCYTE-MONOCYTE.PROGENITOR.CELL","M1.MACROPHAGE.CELL","CMP.CELL",
                   "HEMATOPOIETIC.STEM.CELL","MPP.CELL","HSC.CELL","MONOCYTE.DC.CELL","S100A+")  #"S100A+" è un caso strano

myeloid_notmerged = c("CDC1.CELL","COMMON.LYMPHOID.PROGENITOR.CELL","COMMON.MYELOID.PROGENITOR.CELL","CYCLING.CDC2.CELL",
                      "CYCLING.MYELOID.CELL","DENDRITIC.CELLS.HPA.CELL","LANGERHANS.CELLS.HPA.CELL","MACROPHAGES.HPA.CELL",
                      "MAST.CELL","MEGAKARYOCYTE.CELL","MEGAKARYOCYTE.PROGENITOR.CELL","MHC.CLASS.II.HIGH.MONOCYTE.DC.MACROPHAGE.CELL",
                      "PLACENTA.CAMP+.ANTIMICROBIAL.MONOCYTE.CELL","PLACENTA.DEFENSIN+.ANTIMICROBIAL.MONOCYTE.CELL","PRE-PDC..LYMPHOID.ORIGIN..CELL",
                      "PROLIFERATION.MAST.CELL","PROLIFERATION.MYELOID.CELL","PDC.CELL","ITGAX+.MACROPHAGE.CELL")

## neutrophils_granulosa ##
neutrophils_merged = c("NEUTROPHIL.CELL","NEUTROPHILS.HPA.CELL","GRANULOCYTE.CELL")
neutrophils_notmerged = c("METAMYELOCYTE.BAND.NEUTROPHIL.CELL")

## Alveolar macrophage ##
alv_macro_merged = c("ALVEOLAR.MACROPHAGE.CELL","S.PHASE.MYELOID.CELL","G2.M.PHASE.MYELOID.CELL")

## Tissue-res-macrophage ##
ts_macro_merged = c("LYVE1.MACROPHAGE.CELL","CYCLING.LYVE1.MACROPHAGE.CELL","CDC.CELL","DIAPH3+.MACROPHAGE.CELL")

## pre-pdc myeloid ##
pre_pdc_myeloid = c("MYELOID.PRE-PDC.CELL","PRE-PDC..MYELOID.ORIGIN..CELL")

## hofbauer ##
hofbauer_merged = c("PROLIFERATION.MACROPHAGE.CELL","HOFBAUER.CELLS.HPA.CELL")
hofbauer_notmerged = c("HOFBAUER.CELL")

## BlymphProg ##
blymph_pro_merged = c("PREB.CELL..CYCLING..CELL","PROB.CELL..CYCLING..CELL")
blymph_pro_notmerged = c("LYMPHOID.PRE-PDC.CELL")

## cycling DNDP t cell ##
cycling_DNDP_tcell_merged = c("PROLIFERATION.DN.DP.T.CELL","CYCLING.DP.T.CELL")

## dark ligh GC b cell ##
dark_light_GC_b_cell_merged = c("DARK.ZONE.GC.B.CELL","LIGHT.ZONE.GC.B.CELL")

## breast_milk_macro ##
breast_milk_macro_merged = c("LACTOFERRIN+.MACROPHAGE.CELL","METALLOTHIONEIN+.MACROPHAGE.CELL","INF.RESPONSED.MACROPHAGE.CELL")

## prolif t cell ##
prolif_t_cell_merged = c("PROLIFERATION.T.CELL","TFH.CELL")

## double positive t ##
double_pos_t_merged = c("ABT..ENTRY..CELL","DOUBLE.POSITIVE.T.CELL")

ordered_cells_to_match = c(tcells_merged,bcells_merged,myeloid_merged,neutrophils_merged,alv_macro_merged,
                           ts_macro_merged,pre_pdc_myeloid,hofbauer_merged,blymph_pro_merged,cycling_DNDP_tcell_merged,
                           dark_light_GC_b_cell_merged,breast_milk_macro_merged,prolif_t_cell_merged,double_pos_t_merged,
                           tcells_notmegred,bcells_notmerged,myeloid_notmerged,neutrophils_notmerged,hofbauer_notmerged,
                           blymph_pro_notmerged)

# sum(length(tcells_merged),length(tcells_notmegred),length(bcells_merged),length(bcells_notmerged),length(myeloid_merged),
#     length(myeloid_notmerged),length(neutrophils_merged),length(neutrophils_notmerged),length(alv_macro_merged),length(ts_macro_merged),
#     length(pre_pdc_myeloid),length(hofbauer_merged),length(hofbauer_notmerged),length(blymph_pro_merged),length(blymph_pro_notmerged),
#     length(cycling_DNDP_tcell_merged),length(dark_light_GC_b_cell_merged),length(breast_milk_macro_merged),length(prolif_t_cell_merged),length(double_pos_t_merged))

# length(ordered_cells_to_match)

nnnn = rep(NA,length(ordered_cells_to_match))
for (y in 1:length(nnnn)){
  nnnn[y]=paste(paste("\\<",ordered_cells_to_match[y],sep=""),"\\>",sep="")
}

togrep = paste(nnnn,collapse = "|")
immune_conf = scaled_conf[grepl(togrep,rownames(scaled_conf)),grepl(togrep,colnames(scaled_conf))]
new_names = ordered_cells_to_match[ !ordered_cells_to_match == "S100A+"]
ordered_immune_conf = immune_conf[new_names,new_names]
immune_new_row_names = str_to_title(tolower(rownames(ordered_immune_conf)))
immune_new_col_names = str_to_title(tolower(colnames(ordered_immune_conf)))
rownames(ordered_immune_conf)=immune_new_row_names
colnames(ordered_immune_conf)=immune_new_col_names

adj_immune_names = c("Proliferation immune","Gamma Delta T cell","CD16 Nk cell","NK cell","Cycling T cell",
                     "GZMB CD8 T cell","Cytotoxic CD8 T cell","Cycling T NK cell","CD56 NK cell",
                     "GZMK+ IL7R- CD8 T cell","GZMK NK cell","Tissue Resident NK cell","NK-Cells (HPA)",
                     "T NK cell","Memory T cell","GZMK+ IL7R+ CD8 T cell","KLRB1 CD8 T cell","CD8 T cell",
                     "Memory CD8 T cell","GZMK CD8 T cell","IFN-Activated T cell","Proliferation T NK cell",
                     "CXCL13 Exhausted CD8 T cell","CD4 T cell","Treg cell","KLRB1 Cytotoxic CD4 T cell",
                     "GZMK cytotoxic CD4 T cell", "Effector memory CD4 T cell","Memory CD4 T cell",
                     "Naive CD4 T cell","Naive T cell","INF responsed T cell","Differentiating Treg cell",
                     "Naive CD8 T cell","Agonist T cell","CD8AA+ T cell","GDT cell","Mait cell",
                     "T-Cells (HPA)","Pre-GC B cell","GC committed cell","INF responsed Naive B cell",
                     "Non-Class switch Memory B cell","DZ LZ B cell","Class switch Memory B cell",
                     "Memory B precursor","Naive B cell","LZ B cell","FCRL4 Memory B cell",
                     "S-phase DZ B cell","G2M DZ B cell","Histone high S-phase DZ B cell",
                     "Non-proliferative GC B cell","Plasma cell precursor","PIF1+ G2M DZ B cell","B cell",
                     "Memory B cell","Monocyte cell","CD14 Monocyte","MHCII high CD14 Monocyte",
                     "MHCII low CD14 Monocyte","CD16 Monocyte","CD14 CD16 Monocyte","CDC2 cell",
                     "CD14+ MHCII high Monocyte","CD14+ MHCII low Monocyte","Dendritic cell","Macrophage",
                     "S100A+ Preneutrophil Cycling cell","Myelocyte","Cycling S100A+ Preneutrophil cell",
                     "S100A+ Preneutrophil cell","Promyelocyte","Monocytes (HPA)","CDP cell","GMP cell",
                     "Granulocyte-Monocyte progenitor","M1 Macrophage","CMP cell","Hematopoietic Stem cell",
                     "MPP cell","HSC cell","Monocyte DC cell","Neutrophil","Neutrophils (HPA)","Granulocyte",
                     "Alveolar Macrophage","S-phase Myeloid cell","G2M-phase Myeloid cell","LYVE1 Macrophage",
                     "Cycling LYVE1 Macrophage","CDC cell","DIAPH3+ Macrophage","Myeloid pre-PDC cell",
                     "Pre-PDC Myeloid origin cell","Proliferation Macrophage","Hofbauer cells (HPA)",
                     "PreB cell cycling","ProB cell cycling","Proliferation DN DP T cell",
                     "Cycling DP T cell","Dark zone GC B cell","Light zone GC B cell",
                     "Lactoferrin+ Macrophage","Metallothionein+ Macrophage","INF responsed Macrophage",
                     "Proliferation T cell","Tfh cell","AbT Entry cell","Double positive T cell",
                     "Double negative T cell","T cell","B-Cells (HPA)","PreB cell","ProB cell","CDC1 cell",
                     "Common Lymphoid progenitor","Common Myeloid progenitor","Cycling CDC2 cell",
                     "Cycling Myeloid cell","Dendritic cells (HPA)","Langerhans cells (HPA)",
                     "Macrophages (HPA)","Mast cell","Megakaryocyte","Megakaryocyte progenitor",
                     "MHC class II high Monocyte DC Macrophage","Placenta CAMP+ Antimicrobial Monocyte",
                     "Placenta DEFENSIN+ Antimicrobial Monocyte","Pre-PDC Lymphoid origin cell",
                     "Proliferation Mast cell","Proliferation Myeloid cell","PDC cell","ITGAX+ Macrophage",
                     "Metamyelocyte band neutrophil","Hofbauer cell","Lymphoid pre-PDC cell")    

row.names(ordered_immune_conf)=adj_immune_names
colnames(ordered_immune_conf)=adj_immune_names
## BARS ## "#014636" 
bars = barplot_matrix(TrueLabelsPath,PredLabelsPath)
filt_bars = log10(data.frame("Rel. Aboundande"=bars[new_names,]))
rownames(filt_bars)=str_to_title(new_names)
ha_list = rowAnnotation(" "= anno_barplot(filt_bars, gp = gpar(fill = "#00A087FF"),
                                          bar_width = 1, width = unit(2, "cm")))
# Heatmap # "#a6bddb","#02818a","#016c59"
ordered_immune_conf[ordered_immune_conf==0] <- NA
png("immune_initial_labels_plot.png",width= 10,height= 6,units= "in",res=1200)
PLOT1 = Heatmap(ordered_immune_conf,cluster_columns = F,cluster_rows = F,col=colorRampPalette(c("#E0F2F1","#00796B","#004D40"))(100),
                column_names_rot = 45,na_col = "#fff7fb",row_names_gp = gpar(fontsize = 3,fontface = "bold"),column_names_gp = gpar(fontsize = 3,fontface = "bold"),
                heatmap_legend_param = list(title = "",legend_direction = "horizontal",legend_width = unit(5, "cm")),row_names_side = c("left"),
                right_annotation = ha_list)

draw(PLOT1, heatmap_legend_side = "bottom")
dev.off()

rm(list=setdiff(ls(), c("confusion_matrix","barplot_matrix")))

#### immune system after merging and with Multiassigned, assigned and unassigend categories ####
TrueLabelsPath = "last_merge/try_02/original.csv"
PredLabelsPath = "last_merge/try_02/annotated.csv"
scaled_conf=confusion_matrix(TrueLabelsPath,PredLabelsPath)

tcells_notmegred = c("DOUBLE.NEGATIVE.T.CELL","T.CELL")
bcells_notmerged = c("B-CELLS.HPA.CELL","PREB.CELL","PROB.CELL")
myeloid_notmerged = c("CDC1.CELL","COMMON.LYMPHOID.PROGENITOR.CELL","COMMON.MYELOID.PROGENITOR.CELL","CYCLING.CDC2.CELL"
                      ,"DENDRITIC.CELLS.HPA.CELL","LANGERHANS.CELLS.HPA.CELL","MACROPHAGES.HPA.CELL",
                      "MAST.CELL","MEGAKARYOCYTE.CELL","MEGAKARYOCYTE.PROGENITOR.CELL","MHC.CLASS.II.HIGH.MONOCYTE.DC.MACROPHAGE.CELL",
                      "PLACENTA.CAMP+.ANTIMICROBIAL.MONOCYTE.CELL","PLACENTA.DEFENSIN+.ANTIMICROBIAL.MONOCYTE.CELL","PRE-PDC..LYMPHOID.ORIGIN..CELL",
                      "PROLIFERATION.MAST.CELL","PDC.CELL","ITGAX+.MACROPHAGE.CELL")
neutrophils_notmerged = c("METAMYELOCYTE.BAND.NEUTROPHIL.CELL")
hofbauer_notmerged = c("HOFBAUER.CELL")
blymph_pro_notmerged = c("LYMPHOID.PRE-PDC.CELL")


merged = c("MULTIFUNCTIONAL.T.CELL.2","MULTIFUNCTIONAL.B.CELL","MYELOID.CELL.2","NEUTROPHIL.GRANULOGYTE.CELL",
           "ALVEOLAR.MACROPHAGE.CELL","TISSUE.RESIDENT.MACROPHAGE.CELL","PROLIFERATION.HOFBAUER.CELL",
           "B.LYMPHOCYTE.PROGENITOR.CYCLING.CELL","CYCLING.DN.DP.T.CELL","DARK.LIGHT.ZONE.GC.B.CELL","BREAST.MILK.MACROPHAGE.CELL",
           "DOUBLE.POSITIVE.CELL","PROLIFERATION.DIFFERENTING.T.CELL")#,"CYCLING.MYELOID.CELL","PRE-PDC.MYELOID.CELL","PROLIFERATION.MYELOID.CELL")

ordered_cells_to_match = c(merged,tcells_notmegred,bcells_notmerged,myeloid_notmerged,neutrophils_notmerged,
                           hofbauer_notmerged,blymph_pro_notmerged)

ordered_cells_to_match
nnnn = rep(NA,length(ordered_cells_to_match))
for (y in 1:length(nnnn)){
  nnnn[y]=paste(paste("\\<",ordered_cells_to_match[y],sep=""),"\\>",sep="")
}

togrep = paste(nnnn,collapse = "|")

immune_conf = scaled_conf[grepl(togrep,rownames(scaled_conf)),grepl(togrep,colnames(scaled_conf))]
new_names = ordered_cells_to_match[! ordered_cells_to_match %in% c("CYCLING.MYELOID.CELL", "PRE-PDC.MYELOID.CELL",
                                                                   "PROLIFERATION.MYELOID.CELL","PLACENTA.CAMP+.ANTIMICROBIAL.MONOCYTE.CELL",
                                                                   "PLACENTA.DEFENSIN+.ANTIMICROBIAL.MONOCYTE.CELL","ITGAX+.MACROPHAGE.CELL")]
ordered_immune_conf = immune_conf[new_names,new_names]
immune_new_row_names = str_to_title(tolower(rownames(ordered_immune_conf)))
immune_new_col_names = str_to_title(tolower(colnames(ordered_immune_conf)))
rownames(ordered_immune_conf)=immune_new_row_names
colnames(ordered_immune_conf)=immune_new_col_names

adj_name_pruned_immune = c("Multifunctional T cell","Multifunctional B cell","Myeloid cell",
                           "Neutrophil Granulogyte","Alveolar Macrophage","Tissue Resident Macrophage",
                           "Proliferation Hofbauer cell","B Lymphocyte progenitor cycling cell",
                           "Cycling DN DP T cell","Dark light zone GC B cell","Breast milk Macrophage",
                           "Double positive T cell","Proliferation Differenting T cell",
                           "Double negative T cell","T cell","B-Cells (HPA)","PreB cell","ProB cell",
                           "CDC1 cell","Common Lymphoid progenitor","Common Myeloid progenitor",
                           "Cycling CDC2 cell","Dendritic cells (HPA)","Langerhans cells (HPA)",
                           "Macrophages (HPA)","Mast cell","Megakaryocyte","Megakaryocyte progenitor",
                           "MHC class II high Monocyte DC Macrophage","Pre-PDC Lymphoid origin cell",
                           "Proliferation Mast cell","PDC cell","Metamyelocyte band neutrophil",
                           "Hofbauer cell","Lymphoid pre-PDC cell")
row.names(ordered_immune_conf)=adj_name_pruned_immune
colnames(ordered_immune_conf)=adj_name_pruned_immune

toMatch <- c(new_names,"multiassigned","unassigned")
mmmmm = rep(NA,length(toMatch))
for (yy in 1:length(mmmmm)){
  mmmmm[yy]=paste(paste("\\<",toMatch[yy],sep=""),"\\>",sep="")
}

final_tomatch = paste(mmmmm,collapse = "|")
m4barplot = scaled_conf[grepl(final_tomatch,rownames(scaled_conf)),grepl(final_tomatch,colnames(scaled_conf))]
bars = data.frame("Assigned"=1-apply(m4barplot[,c("multiassigned","unassigned")],1,sum),
                  "Multiassigned"= m4barplot["multiassigned"],"Unassigned"=m4barplot["unassigned"])
bars = bars[new_names,]                                                                                                                                                                      

ha_list = rowAnnotation(" "= anno_barplot(bars, gp = gpar(fill = c("#fdbf6f","#cab2d6","#b2df8a")), 
                                          bar_width = 1, width = unit(2, "cm")))

png("cleaned_lists_merged_immune_plot.png",width= 10,height= 6,units= "in",res=1200)
PLOT2 <- Heatmap(ordered_immune_conf,cluster_columns = F,cluster_rows = F,show_row_dend = FALSE,
                 show_column_dend = FALSE,row_names_gp = gpar(fontsize = 8,fontface = "bold"),column_names_gp = gpar(fontsize = 8,fontface = "bold"),column_names_side = c("bottom"),
                 row_names_side = c("left"),right_annotation = ha_list,column_names_rot = 45,na_col = "#fff7fb",
                 heatmap_legend_param = list(title = ""),col=colorRampPalette(c("#E0F2F1","#00796B","#004D40"))(100))#col = colorRamp2(c(0, 0.5, 1),c("#f7f7f7","#c2a5cf","#762a83"))
lgd = Legend(labels = c("Assigned","Multiassigned","Unassigned"), title = "Class", legend_gp = gpar(fill = c("#fdbf6f","#cab2d6","#b2df8a")))
draw(PLOT2,annotation_legend_list = lgd) 
dev.off()

rm(list=setdiff(ls(), c("confusion_matrix","barplot_matrix")))

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
}
