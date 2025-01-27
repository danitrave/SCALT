library(ggplot2)
library(viridis)
original = read.table("68k_pbmc_barcodes_annotation.tsv",sep="\t",row.names = NULL,header = 1)
colnames(original)<-c("TSNE.1","TSNE.2","barcodes","Types")

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

p_original <-ggplot(original, aes(x=TSNE.1, y=TSNE.2,color=Types)) + 
  geom_point(size=1)+scale_colour_manual(values = c25)+#scale_color_brewer(palette = "Spectral")+
  theme(legend.key.size=unit(5,"point"))+ guides(colour = guide_legend(override.aes = list(size=5)))

ggsave("zhang68k_originalTSNE.png",device = "png", width = 30, height = 20, units = "cm")

