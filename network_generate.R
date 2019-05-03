## Load edge weight and node weight file
library(igraph)
edge.weight = read.table("C:/Users/fyan/EW_dmGWAS/1_29/edge_weight_scale.txt", header=F, as.is=T)
node.weight = read.table("C:/Users/fyan/EW_dmGWAS/1_29/Asian_node_weight_scale.txt", as.is=T)

## Generate graph
G = graph.data.frame(edge.weight[,1:2], directed=F)
E(G)$weight = edge.weight[, 3]
z_node = node.weight[,2]
names(z_node) = node.weight[,1]
V(G)$weight = z_node[V(G)$name]

## Load output file generated from EW_dmGWAS
res.mat = read.delim("C:/Users/fyan/EW_dmGWAS/1_29/Asian.lambda_1_permutated.txt", as.is=T)
res.mat = res.mat[order(res.mat[, "z_perm"], decreasing=T), ]
head(res.mat)

## Top 1 module
top1_genes = strsplit(res.mat[1,1], split=" ")[[1]]
top1_genes
top1_subG <- induced.subgraph(G, top1_genes)
vcount(top1_subG)
top1_edge = cbind(get.edgelist(top1_subG), E(top1_subG)$weight)
top1_node = cbind(V(top1_subG)$name, V(top1_subG)$weight)
write.table(top1_edge,file="top1_edge.txt", row.names=FALSE, sep="\t")
write.table(top1_node,file="top1_node.txt", row.names=FALSE, sep="\t")

## All significant modules
sig.res.mat <- res.mat[res.mat$z_perm > 1.96, ]
module.genes <- sapply(sig.res.mat[,1], function(u)strsplit(u, split=" ")[[1]] )
module.genes <- unique(unlist(module.genes))
length(module.genes)
All_subG <- induced.subgraph(G, All_genes)
vcount(All_subG)
All_edge = cbind(get.edgelist(All_subG), E(All_subG)$weight)
All_node = cbind(V(All_subG)$name, V(All_subG)$weight)
write.table(All_edge,file="All_edge.txt", row.names=FALSE, sep="\t")
write.table(All_node,file="All_node.txt", row.names=FALSE, sep="\t")
