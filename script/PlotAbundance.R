#!
args <- commandArgs(TRUE)

infile  <- args[1]
outfile <- args[2]

library(ggplot2)

abd <- read.table(infile, header=TRUE, check.names = TRUE)
class(abd$species)

###draw
jpeg(outfile,width=2400,height=2000,res=300)
ggplot(abd, aes(x=species, y=abundance, group = species)) + 
  theme(axis.text.x=element_text(size=7, angle=45)) +
  geom_violin(trim=F,draw_quantiles = c(0.25, 0.5, 0.75),col="grey90") +
  geom_jitter(width=0.2,cex=0.5,col="#db650e")
dev.off()
