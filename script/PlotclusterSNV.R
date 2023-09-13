#!
args <- commandArgs(TRUE)

BAFfile     <- args[1]
ABDfile     <- args[2]
clusterfile <- args[3]
outfile     <- args[4]

###env
require(clusterSim)
library(vegan)
library(proxy)

###function: x is a distance matrix and k the number of clusters
pam.clustering <- function(x, k) {
	require(cluster)
	cluster <- as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
	return(cluster)
}

###input
BAFmat <- read.table(BAFfile, header=TRUE, sep="\t", row.names = 1)
ABDmat <- read.table(ABDfile, header=TRUE, sep="\t", row.names = 1)

###distance
if(ncol(BAFmat) < 3){
	distmat <- dist(BAFmat, method = "euclidean")
}else{
	distmat <- as.dist((1 - cor(t(BAFmat)))/2)
}

###clustering
ROUNDS <- 6
score.cluster <- 0
num.cluster <- NULL
score.max <- 0
for(k in 1:ROUNDS){
	if(k==1){ 
		num.cluster=1
	} else{
		print(k)
		cluster_temp <- pam.clustering(distmat, k)
		score.cluster[k] <- abs(index.G1(BAFmat, cluster_temp, d=distmat, centrotypes="medoids"))
		if(score.cluster[k] > score.max){
			score.max <- score.cluster[k]
			num.cluster <- k
		}
	}
}

###decide the number of cluster
#pdf(clusterfile, width=8, height=8)
jpeg(clusterfile,width=2400,height=2000,res=300)
plot(score.cluster, type="h", xlab="k clusters", ylab="CH Index", main="Optimal number of clusters")
dev.off()

###new BAFmat
cluster_final <- pam.clustering(distmat, k=num.cluster) 
BAFmat.new <-  cbind(BAFmat, cluster_final)

###draw
number <- ncol(BAFmat.new) - 1
#pdf(pdffile,width=10,height=8)
jpeg(outfile,width=2400,height=2000,res=300)

#setting
cols <- c("orange","deepskyblue","#db5353","red","green","purple")
layout(matrix(c(1,2),2,1),heights = c(1:5))
 
#panel A 
par(mar = c(0, 5, 2, 2))
plot("", xlim=c(0,number+1), ylim=c(0,0.1), xlab=NA, ylab="relative abundance", xaxt="n")
for(i in 1:ncol(ABDmat)){
	rect(i-1,0,i, ABDmat[1,i],col="#f4bb82",border="orange")
}

#panel B 
par(mar = c(5, 5, 0, 2)) 
plot("",xlim=c(0,number+1),ylim=c(0,1),xlab="timepoints",ylab="B Allele Frequency",xaxt="n")
for(i in 1:nrow(BAFmat.new)){
	col_draw = cols[BAFmat.new[i, number + 1]]
	x=seq(1,number,by=1)
	lines(x, as.numeric(BAFmat.new[i,1:number]),type="o", lwd=0.1, cex=0.1, col=col_draw)
}
axis(1,at=seq(1,number,by=1),labels=colnames(BAFmat)[1:number])

dev.off()
