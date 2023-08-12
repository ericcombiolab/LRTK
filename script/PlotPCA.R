#！
args <- commandArgs(TRUE)

infile   <- args[1]
metafile <- args[2]
outfile  <- args[3]

library(factoextra)
library(ggforce)
library(vegan)

datgroup <- read.table(metafile,header=TRUE,sep="\t",row.names=1)
group <- datgroup$Group

dat <- read.table(infile,header=TRUE,row.names=1,check.names = FALSE)
datdist <- vegdist(t(dat),method = "bray")

jpeg(outfile,width=2400,height=2000,res=300)
res.pca <- prcomp(datdist,  scale = TRUE)
fviz_pca_ind(res.pca, label='none',repel = TRUE, col.ind=datgroup$Group) +
  ggforce::geom_mark_ellipse(aes(fill = datgroup$Group,color = datgroup$Group)) + theme(legend.position = 'bottom') + coord_equal()
dev.off()
