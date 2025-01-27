library(ggplot2)

probs = read.table("GENERAL_PROBABILITIES.tsv")
entro = read.table("genes_entropy.tsv")
all(row.names(probs)==row.names(entro))
ct_probs = read.table("CELL_TYPES_PROBABILITIES.tsv")
keep = apply(ct_probs,1,max)>0.10
genes2keep = row.names(probs)[keep]

scatter.smooth(probs[genes2keep,"RATIO"],entro[genes2keep,"entropies"],col="red")
cor.test(probs[genes2keep,"RATIO"],entro[genes2keep,"entropies"])
df = data.frame(p=probs[genes2keep,"RATIO"],e=entro[genes2keep,"entropies"])

row.names(df)=genes2keep
upper_bound = quantile(df$e,0.95)
lower_bound = quantile(df$e,0.05)
# Basic scatter plot
PLOT=ggplot(df, aes(x=p, y=e)) + geom_point(color="#DC0000B2") +geom_hline(yintercept = lower_bound)+
  geom_hline(yintercept = upper_bound)+labs(y= "Entropies", x = "Probabilites")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("Figure_5.png",PLOT,dpi = 1200)
dev.off()
dim(df[df$e>lower_bound & df$e<upper_bound,])

