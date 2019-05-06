## Change directory to your directory
setwd("path to your directory")

## Load edge weight file, node weight file, EW_dmGWAS output file
edge.weight = read.table("path to your edge weight file", header=F, as.is=T)
node.weight = read.table("path to your node weight file", as.is=T)
res.mat = read.delim("path to your EW_dmGWAS output file", as.is=T)

res.mat = res.mat[order(res.mat[, "z_perm"], decreasing=T), ]
## All significant modules
sig.res.mat <- res.mat[res.mat$z_perm > 1.96, ]

## Generate graph
library(igraph)
G = graph.data.frame(edge.weight[,1:2], directed=F)
E(G)$weight = edge.weight[, 3]
z_node = node.weight[,2]
names(z_node) = node.weight[,1]
V(G)$weight = z_node[V(G)$name]

## Generate top n modules 
top_modules <- function(x){
top_genes <- sapply(sig.res.mat[1:x,1], function(u)strsplit(u, split=" ")[[1]])
module_genes <- unique(unlist(top_genes))
subG <- induced.subgraph(G, module_genes)
edge <- cbind(get.edgelist(subG), E(subG)$weight)
node <- cbind(V(subG)$name, V(subG)$weight)
write.table(edge,file="edge_top.txt", row.names=FALSE, quote=F, sep="\t")
write.table(node,file="node_top.txt", row.names=FALSE, quote=F, sep="\t")}

## Top 10 modules
top_modules(10)
                    
## Top 50 modules
top_modules(50)                   
