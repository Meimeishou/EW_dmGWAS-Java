#Astrid M Manuel
#03/08/2019

#This script produces edge and node files of top networks resulting from EW_dmGWAS
#files that result are input for graphing top modules in Cytoscape

#read in results of EW_dmGWAS permutation, as well as edge and node files that were input for EW_dmGWAS
setwd("C://Users//amanuel1//Dev//ZhaoLab//MS//GeneExpressionProfiles//EW//PLChronic_and_2011GWAS")
edge.weight <- read.delim("PL_Chronic_Active_EW_shared_2011.txt",header=F,as.is=T)
node.weight <- read.delim("PL_Chronic_Active_NW_shared_2011.txt",header=F,as.is=T)
res.mat <- read.delim("PL_Chronic_Active_vs_Control_2011_permutation.txt",as.is=T)
res.mat = res.mat[order(res.mat[, "z_perm"], decreasing=T), ]
sig.res.mat = res.mat[res.mat$z_perm > 1.96, ]
sapply(sig.res.mat[,1], function(u)strsplit(u, split=" ")[[1]] ) -> module.genes
write.table(head(res.mat, 10), "Top10ModulesSummary.txt", quote =F, row.names=F, sep = "\t")

#number of top modules to map
n <- 10

for(i in seq(1:n)){
  
  #genes of top module #i
  top.genes <- unlist(module.genes[i], use.names = F)
  
  #edges and edges weights for module i
  idx.ew <- which(edge.weight[,1] %in% top.genes & edge.weight[,2] %in% top.genes)
  ew <- edge.weight[idx.ew,]
  ew[,4] <- 2*pnorm(-abs(ew[,3]))
  ew[,5] <- -log10(ew[,4])
  ew.file <- ew[,c(1,2,5)]
  colnames(ew.file) <- c("source", "target", "edge weight")
  write.table(ew.file, paste("EW_PL_Chronic_Active_2011_TopModule", i, ".txt", sep = ""), quote=F, row.names = F, sep = "\t")
  
  #nodes and node weights for module i
  idx.nw <- which(node.weight[,1] %in% top.genes)
  nw <- node.weight[idx.nw,]
  nw[,3] <- 2*pnorm(-abs(nw[,2]))
  nw[,4] <- -log10(nw[,3])
  nw.file <- nw[,c(1,4)]
  colnames(nw.file) <- c("key", "node weight")
  write.table(nw.file, paste("NW_PL_Chronic_Active_2011_TopModule", i, ".txt", sep = ""), quote=F, row.names = F, sep = "\t")

}

#to graph top 5 modules, merged
#gene list of top 5 modules
top.5M.genes <- unlist(module.genes[1:10], use.names = F)
idx.ew <- which(edge.weight[,1] %in% top.5M.genes & edge.weight[,2] %in% top.5M.genes)
ew <- edge.weight[idx.ew,]
ew[,4] <- 2*pnorm(-abs(ew[,3]))
ew[,5] <- -log10(ew[,4])
ew.file <- ew[,c(1,2,5)]
colnames(ew.file) <- c("source", "target", "edge weight")
ew.file <- ew.file[which(ew.file[,3]>1.3), ]
write.table(ew.file, "EW_PL_Chronic_Active_2011_TopModule5M.txt", quote=F, row.names = F, sep = "\t")

top.5M.NW <- unique(c(ew.file[,1], ew.file[,2]))
idx.nw <- which(node.weight[,1] %in% top.5M.NW)
nw <- node.weight[idx.nw,]
nw[,3] <- 2*pnorm(-abs(nw[,2]))
nw[,4] <- -log10(nw[,3])
nw.file <- nw[,c(1,4)]
write.table(nw.file, "NW_PL_Chronic_Active_2011_TopModule5M.txt", quote=F, row.names = F, sep = "\t")

#gene list of all modules
top.genes <- unlist(module.genes[1:10], use.names = F)
idx.ew <- which(edge.weight[,1] %in% top.genes & edge.weight[,2] %in% top.genes)
ew <- edge.weight[idx.ew,]
ew[,4] <- 2*pnorm(-abs(ew[,3]))
ew[,5] <- -log10(ew[,4])
ew.file <- ew[,c(1,2,5)]
colnames(ew.file) <- c("source", "target", "edge weight")
#ew.file <- ew.file[which(ew.file[,3]>1.3), ]
write.table(ew.file, "EW_PL_Chronic_Active_2011_10M.txt", quote=F, row.names = F, sep = "\t")

top.NW <- unique(c(ew.file[,1], ew.file[,2]))
idx.nw <- which(node.weight[,1] %in% top.NW)
nw <- node.weight[idx.nw,]
nw[,3] <- 2*pnorm(-abs(nw[,2]))
nw[,4] <- -log10(nw[,3])
nw.file <- nw[,c(1,4)]
write.table(nw.file, "NW_PL_Chronic_Active_2011_10M.txt", quote=F, row.names = F, sep = "\t")

