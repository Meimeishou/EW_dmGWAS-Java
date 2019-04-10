load("C:\\Users\\fyan\\Desktop\\18cleft\\edge_weight.RData")
ls()
dim(df)
head(df)
df <- na.omit(df)
head(df)
dim(df)
summary(df$e)
scale(df$e) -> df$e_scale
head(df)
write.table(df,"edge_weight.txt",sep="\t",row.names=FALSE)
write.table(df[,c(1,2,9)],"edge_weight_scale.txt",sep="\t",row.names=FALSE)
ls(0
)
ls()
gwas_asian <- read.delim("C:/Users/fyan/Desktop/18cleft/dmGWAS_R/PASCAL_original_10_26/adjusted/GWAS_Asian_adjusted.txt",as.is=T)
head(gwas_asian)
gwas_asian$e <- qnorm(1- gwas_asian$pvalue/2)
head(gwas_asian)
hist(gwas_asian$e)
gwas_asian$e_scale <- scale(gwas_asian$e)
head(gwas_asian)
write.table(gwas_asian,"node_weight.txt",sep="\t",row.names=FALSE)
write.table(gwas_asian[,c(1,4)],"node_weight_scale.txt",sep="\t",row.names=FALSE)
q()
getwd()
library(igraph)
edge.weight = read.table("C:/Users/fyan/EW_dmGWAS/1_29/edge_weight_scale.txt", header=T, as.is=T)
head(edge.weight)
edge.weight = read.table("C:/Users/fyan/EW_dmGWAS/1_29/edge_weight_scale.txt", header=F, as.is=T)
dim(edge.weight)
G = graph.data.frame(edge.weight[,1:2], weight=edge.weight[,3])
G = graph.data.frame(edge.weight[,1:2], directed=F)
E(G)$weight = edge.weight[, 3]
res.mat = read.delim("C:/Users/fyan/EW_dmGWAS/1_29/Asian.lambda_1_permutated.txt", as.is=T)
head(res.mat)
res.mat = res.mat[order(res.mat[, "z_perm"]), ]
head(res.mat)
res.mat = res.mat[order(res.mat[, "z_perm"], decreasing=T), ]
head(res.mat)
genes = strsplit(res.mat[1,1], split=" ")[[1]]
genes
induced.subgraph(G, genes) -> subG
vcount(subG)
cbind(get.edgelist(subG), E(subG)$weight)
cbind(V(subG)$name, V(subG)$weight)
node.weight = read.table("C:/Users/fyan/EW_dmGWAS/1_29/Asian_node_weight_scale.txt", as.is=T)
head(node.weight)
z_node = node.weight[,2]
names(z_node) = node.weight[m1]
names(z_node) = node.weight[,1]
V(G)$weight = z_node[V(G)$name]
 induced.subgraph(G, genes) -> subG
vcount(G)
ecount(G)
summary(E(G)$weight)
summary(V(G)$weight)
dim(node.weight)
length(intersect(node.weight[,1], V(G)$name))
genes
ecount(subG)
get.edgelist(subG)
cbind(get.edgelist(subG), E(subG)$weight)
cbind(V(subG)$name, V(subG)$weight)
sig.res.mat = res.mat[res.mat$z_perm > 1.96, ]
sapply(sig.res.mat[,1], function(u)strsplit(u, split=" ")[[1]] ) -> module.genes
length(module.genes)
dim(res.mat)
dim(sig.res.mat)
module.genes = unique(unlist(module.genes))
length(module.genes)
sort(module.genes)
history()
history()
history(1000)

