#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("pheatmap")
library(RColorBrewer)
library(stringr) 
library("ComplexHeatmap")
library("ggsci")

mypal <- pal_npg("nrc", alpha = 1.0)(9)
#scales::show_col(mypal)
# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF"

genes_heatmap_plot <- function(t,s,o){
  table = read.table(t,sep="\t",header = T,row.names=1)
  original = read.table(o,sep="\t",header = T)
  scalt = read.table(s,sep="\t",header = T)
  cell_types_order = unique(original$CELL_ANNOTATION)
  cell_types_order = gsub("_", ".", cell_types_order)
  
  lists = list.files("../cell_types/")
  
  L = list()
  for (e in lists){
    file2read = paste("../cell_types/",e,sep="")
    l = read.table(file2read,header = F)
    genes_vector = l[,"V1"]
    new_name = gsub("_",".",strsplit(e,"_most",fixed=T)[[1]][1])
    L[[new_name]]=genes_vector
  }
  
  selected_L = L[cell_types_order]
  genes_for_plot_list = list()
  for (j in names(selected_L)){
    A = selected_L[[j]]
    B = Reduce(union,selected_L[!names(selected_L) == j])
    I = intersect(A,B)
    U = A[!A %in% I]
    genes_for_plot_list[[j]]=U
  }
  
  bar = scalt$CELL_ANNOTATION
  uniq_annot = unique(bar)
  color_list = list()
  for (w in uniq_annot){
    COLOR = "#"
    if (w == "multiassigned"){
      color_list[["multiassigned"]]="#f0f0f0"
    }else{
      for (r in 1:6){
        rd = sample(c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"),1,replace = TRUE)
        COLOR = paste(COLOR,rd,sep="")
      }
      color_list[[w]]=COLOR
    }
  }
  
  final_colors = list("Cell types"=unlist(color_list))
  ha = HeatmapAnnotation("Cell types" = bar,col = final_colors,gp = gpar(fontsize = 1),
                         annotation_legend_param = list(labels_gp=gpar(fontsize = 5)))
  LLL = names(color_list)[! names(color_list) == "multiassigned"]
  split = rep(LLL,as.vector(lengths(genes_for_plot_list)))
  GENES = unlist(genes_for_plot_list,use.names = F)
  filtered_table = table[GENES,]
  filtered_table[filtered_table >= 1] <- 1 #  #0570b0  
  #"#fff5eb","#3690c0"
  PLOT1 = Heatmap(filtered_table,cluster_columns = F,cluster_rows = F,
                  use_raster=F, col = col = c("white","#C51162"),show_column_names = F,
                  show_row_names = F,top_annotation = ha,show_heatmap_legend = F,
                  row_split = split,row_title = NULL,border_gp = gpar(col = "#607D8B"),border=T)
  return(PLOT1)
}

#### Neurons ####
t = "../neuron_data_for_plot/neuron_table_complete.tsv"
s = "../neuron_data_for_plot/neuron_scalt_anno_acomplete.tsv"
o = "../neuron_data_for_plot/neuron_original_anno_complete.tsv"

neuron_plot = genes_heatmap_plot(t,s,o)
output_name = "genes_heatmap_plot_neurons.png"
png(output_name,width= 12,height= 8,units= "in",res=1200)
draw(neuron_plot,heatmap_legend_side = "bottom")
dev.off()

#### Fibroblasts ####
t = "../fibroblast_data_for_plot/fibroblast_table_complete.tsv"
s = "../fibroblast_data_for_plot/fibroblast_scalt_anno_acomplete.tsv"
o = "../fibroblast_data_for_plot/fibroblast_original_anno_complete.tsv"

fibro_plot = genes_heatmap_plot(t,s,o)
output_name = "genes_heatmap_plot_fibroblasts.png"
png(output_name,width= 12,height= 8,units= "in",res=1200)
draw(fibro_plot,heatmap_legend_side = "bottom")
dev.off()

##### Immune system ####
t = "../immune_data_for_plot/immune_table_complete.tsv"
s = "../immune_data_for_plot/immune_scalt_anno_acomplete.tsv"
o = "../immune_data_for_plot/immune_original_anno_complete.tsv"

immune_plot = genes_heatmap_plot(t,s,o)
output_name = "genes_heatmap_plot_immune.png"
png(output_name,width= 12,height= 8,units= "in",res=1200)
draw(immune_plot,heatmap_legend_side = "bottom")
dev.off()

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
}