#obtaining edge weights for PL Chronic Active MS vs Control

# control PPI values
Control_GEmatrix <- nodupGeneRecords[,1:10]
Control_coexpression <- cor(t(Control_GEmatrix))
network <- read.table("C:\\Users\\amanuel1\\Dev\\ZhaoLab\\MS\\HPRD_Release9_062910\\network.txt", as.is=T)
Control_PPI_pccValues <- data.frame()
Control_genes.with.expression = rownames(Control_coexpression)

for(i in seq(1:nrow(network))){
  
  is.element(network[i, 1], genes.with.expression) -> check1
  is.element(network[i, 2], genes.with.expression) -> check2
  
  if( !check1 | (!check2) )next
  
  Control_PPI_pccValues[i,1] <- network[i, 1]
  Control_PPI_pccValues[i,2] <- network[i, 2]
  Control_PPI_pccValues[i,3] <- Control_coexpression[network[i, 1],network[i, 2]]
}
Control_PPI_pccValues_omitNA <- na.omit(Control_PPI_pccValues)

Control_PPI_pccValues_rmNA_rmSelfInteractions <- Control_PPI_pccValues_omitNA[which(Control_PPI_pccValues_omitNA[,1]!=Control_PPI_pccValues_omitNA[,2]), ]


# PL chronic active PPI values
PL_ChronicActive_GEmatrix <- nodupGeneRecords[,18:24]
PL_ChronicActive_coexpression <- cor(t(PL_ChronicActive_GEmatrix))
PL_ChronicActive_PPI_pccValues <- data.frame()
PL_ChronicActive_genes.with.expression = rownames(PL_ChronicActive_coexpression)

for(i in seq(1:nrow(network))){
  
  is.element(network[i, 1], PL_ChronicActive_genes.with.expression) -> check1
  is.element(network[i, 2], PL_ChronicActive_genes.with.expression) -> check2
  
  if( !check1 | (!check2) )next
  
  PL_ChronicActive_PPI_pccValues[i,1] <- network[i, 1]
  PL_ChronicActive_PPI_pccValues[i,2] <- network[i, 2]
  PL_ChronicActive_PPI_pccValues[i,3] <- PL_ChronicActive_coexpression[network[i, 1],network[i, 2]]
}
PL_ChronicActive_PPI_pccValues_omitNA <- na.omit(PL_ChronicActive_PPI_pccValues)

PL_ChronicActive_PPI_pccValues_rmNA_rmSelfInteractions <- PL_ChronicActive_PPI_pccValues_omitNA[which(PL_ChronicActive_PPI_pccValues_omitNA[,1]!=PL_ChronicActive_PPI_pccValues_omitNA[,2]), ]

# matching case vs control shared edges for edge-weight calculations

control_PCC_values <- Control_PPI_pccValues_omitNA_rmSelfInteractions
case_PCC_values <- PL_ChronicActive_PPI_pccValues_rmNA_rmSelfInteractions

Ncontrol <- 10
Ncase <- 7

apply(control_PCC_values, 1, function(u){
  paste( sort(c(u[1], u[2])), collapse="|")
}) -> ctrl.str

apply(case_PCC_values, 1, function(u){
  paste(sort(c(u[1], u[2])), collapse="|")
}) -> case.str

shared.edge.str = intersect(case.str, ctrl.str)
match(shared.edge.str, case.str) -> idx.1
case_PCC_values.shared = case_PCC_values[idx.1, ]

match(shared.edge.str, ctrl.str) -> idx.2
control_PCC_values.shared = control_PCC_values[idx.2, ]


# Obtaining F(x): Fisher's transformation (Jia et al., 2014)
# Built in log() performs natural logarithm - NOT log base 10 as one would expect
F.equation = function(x){
  (1/2) * log((1+x)/(1-x))
}

# Calculating edge-weights
X <- data.frame()
denominator = sqrt((1/(Ncase-3)) + (1/(Ncontrol-3))) 
for(k in 1:nrow(case_PCC_values.shared)){
  X[k,1] <- case_PCC_values.shared[k,1]
  X[k,2] <- case_PCC_values.shared[k,2]
  r.case <- case_PCC_values.shared[k,3]
  r.control <- control_PCC_values.shared[k,3]
  f.case = F.equation(r.case)
  
  f.control = F.equation(r.control)
  
  X[k,3] = (f.case - f.control)/denominator
  X[k,4] = f.case
  X[k,5] = f.control
}

EW <- scale(X[,3])
EW <- abs(EW)
EW <- qnorm(1-(pnorm(EW,lower.tail=F)*2))
EW_file <- cbind(X[,1],X[,2],EW)
write.table(EW_file, "PL_Chronic_Active_EW.txt", row.names = F, col.names = F, quote = F)

#getting the node weights and matching to ew
node_weights <- read.table("C:\\Users\\amanuel1\\Dev\\ZhaoLab\\MS\\GeneExpressionProfiles\\EW\\2011_original_NW.txt", as.is = T)
shared.genes <- intersect(c(EW_file[,1],EW_file[,2]),node_weights[,1])
shared.nodes.weight <- node_weights[match(shared.genes, node_weights[,1]), ]
qnorm(shared.nodes.weight[,2]/2, lower.tail = F) -> node.z
shared.nodes.weight$z <- node.z

NW_file <- cbind(shared.nodes.weight[,1], shared.nodes.weight[,3])
write.table(EW_file, "PL_Chronic_Active_NW.txt", row.names = F, col.names = F, quote = F)


