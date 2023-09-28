#!
args <- commandArgs(TRUE)

snvfile  <- args[1]
dbfile   <- args[2]
outfile  <- args[3]
barcode1 <- args[4]
barcode2 <- args[5]

###input
snps <-read.table(snvfile,header=TRUE)
genome <- read.table(dbfile)

###genomic information
chr.lens <- genome[,1]
chr.w <- chr.lens / sum(chr.lens)
chr.offsets <- c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1) 

###function
SetUpPlot <- function(y.lab, y.min, y.max, x.lab, lab.chr){
  chr.lens <- genome[,1]
  chr.w <- chr.lens / sum(chr.lens)
  suppressWarnings( par( mar=c(1.5,3.6,0.1,0),mgp=c(2,-0.2,-1.1) ) ) 
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(y.min, y.max), xlab = "", ylab = "", main = "", xaxt = "n")
  mtext(side=1, line=0.5, x.lab, cex=0.75, outer=FALSE)
  mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)
  ww <- as.vector(rbind(chr.w, chr.w)) / 2
  lab.vals <- (c(1:length(chr.w)))
  odd.ix <- lab.vals %% 2 == 1
 
  chr.offsets <- c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1)
  for (i in 1:(length(chr.offsets) - 1)) {
    use.col <- ifelse(i%%2 == 1, "grey90", "grey98")
    rect(xleft = chr.offsets[i], ybottom = y.min, xright = chr.offsets[i + 1], ytop = y.max, col = use.col, border = NA)
  }
  chr.dat = list(chr.w, chr.lens, chr.offsets)
}

#pdf(outfile,width=12,height=6)
jpeg(outfile,width=2400,height=2000,res=300)
##setting
snp.cola  <- "deepskyblue"
snp.colb  <- "orange"
par(mfrow=c(2,1), cex=0.75, las=1)
par(oma=c(0.1,2,0.1,2)) 

SetUpPlot("BAF", 0, 1, barcode1, T)
for(i in 1:nrow(snps)){
	a <- strsplit(as.character(snps[i,1]), "_")
	b <- strsplit(a[[1]][3], ":")
	chr <- as.numeric(b[[1]][1])
	snp.crds  <- as.numeric(b[[1]][2])
	snp.aaf   <- 1-snps[i,2]
	snp.baf   <- snps[i,2]
	genome.crds <- chr.offsets[chr] + snp.crds / chr.lens[chr] * chr.w[chr]
	points(genome.crds, snp.aaf, col=snp.cola, pch=16, cex=0.8)
	points(genome.crds, snp.baf, col=snp.colb, pch=16, cex=0.8)
}

SetUpPlot("BAF", 0, 1, barcode2, T)
for(i in 1:nrow(snps)){
	a <- strsplit(as.character(snps[i,1]), "_")
	b <- strsplit(a[[1]][3], ":")
	chr <- as.numeric(b[[1]][1])
	snp.crds  <- as.numeric(b[[1]][2])
	snp.aaf   <- 1-snps[i,3]
	snp.baf   <- snps[i,3]
	genome.crds <- chr.offsets[chr] + snp.crds / chr.lens[chr] * chr.w[chr]
	points(genome.crds, snp.aaf, col=snp.cola, pch=16, cex=0.8)
	points(genome.crds, snp.baf, col=snp.colb, pch=16, cex=0.8)
}

dev.off()
