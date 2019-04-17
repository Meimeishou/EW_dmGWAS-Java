
#Calculating edge weights for EW_dmGWAS input, matching PPI edges to node weights.
#04/10/2019

#gene-level p-values
node_weights <- read.table("C:/Users/amanuel1/Dev/ZhaoLab/MS/GeneExpressionProfiles/EW/MS_EUR_2018/NW_MS_EUR_2018.txt", as.is = T)

#numeric matrix of gene expression values for control samples (gene IDs as rownames):
GE.control <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CONTROL", as.is = T)

#numeric matrix of gene expression values for control samples (gene IDs as rownames):
GE.case <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CASE", as.is = T)

#Protein-prontein interaction network (2 columns for each interaction between 2 genes)
PPI <- read.table("PATH_TO_PPI", as.is = T)


#Obtaining PCC
obtainPCC <- function(GE.matrix){
  PCCvalues <- cor(t(GE.matrix))
  return(PCCvalues)
}

#matching to PPI
matchPCCtoPPI <- function(PCCvalues, PPI){
  
  genes.with.expression <- row.names(GE.case)
  
  PPI_PCC <- data.frame()
  
  
  for(i in seq(1:nrow(PPI))){
    
    is.element(PPI[i, 1], genes.with.expression) -> check1
    is.element(PPI[i, 2], genes.with.expression) -> check2
    
    if( !check1 | (!check2) )next
    
    PPI_PCC[i,1] <- PPI[i, 1]
    PPI_PCC[i,2] <- PPI[i, 2]
    PPI_PCC[i,3] <- PCCvalues[PPI[i, 1],PPI[i, 2]]
  }
  
  PPI_PCC_rmNA <- na.omit(PPI_PCC)
  
  PPI_PCC_rmNA_rmSI <- PPI_PCC_rmNA[which(PPI_PCC_rmNA[,1]!=PPI_PCC_rmNA[,2]), ]
  
  return(PPI_PCC_rmNA_rmSI)
}


# Obtaining F(x): Fisher's transformation (Jia et al., 2014)
F.equation = function(x){
  (1/2) * log((1+x)/(1-x))
}



message("Obtaining CASE co-expression values...")
casePCC <- obtainPCC(GE.case)
message("Obtaining CONTROL co-expression values...")
controlPCC <- obtainPCC(GE.control)

message("Iterating through PPI network and matching co-expression values...")
control_PPI_PCC <- matchPCCtoPPI(controlPCC, PPI)
case_PPI_PCC <- matchPCCtoPPI(casePCC, PPI)  

# matching case vs control shared edges for edge-weight calculations
message("Matching case and control shared edges...")
apply(control_PPI_PCC, 1, function(u){
    paste( sort(c(u[1], u[2])), collapse="|")
  }) -> ctrl.str

  apply(case_PPI_PCC, 1, function(u){
    paste(sort(c(u[1], u[2])), collapse="|")
  }) -> case.str

  shared.edge.str = intersect(case.str, ctrl.str)
  
  match(shared.edge.str, ctrl.str) -> idx.1
  control_PPI_PCC.shared = control_PPI_PCC[idx.1, ] 
  
  match(shared.edge.str, case.str) -> idx.2
  case_PPI_PCC.shared = case_PPI_PCC[idx.2, ]
  
  
  
# obtaining edge weights
  
message("Calculating edge weights...")
  Ncase <- ncol(GE.case)
  Ncontrol <- ncol(GE.control)
  
  denominator = sqrt((1/(Ncase-3)) + (1/(Ncontrol-3))) 
  
  X <- data.frame()
  
  for(k in 1:nrow(case_PPI_PCC.shared)){
    X[k,1] <- case_PPI_PCC.shared[k,1]
    X[k,2] <- case_PPI_PCC.shared[k,2]
    r.case <- case_PPI_PCC.shared[k,3]
    r.control <- control_PPI_PCC.shared[k,3]
    f.case = F.equation(r.case)
    
    f.control = F.equation(r.control)
    
    X[k,3] = (f.case - f.control)/denominator
    X[k,4] = f.case
    X[k,5] = f.control
  }

X <- X[-c(which(is.infinite(X[,3]))),]
EW <- X[,3]
EW <- scale(EW)
EW <- abs(EW)
EW <- qnorm(1-(pnorm(EW,lower.tail=F)*2))
X$EW <- EW
EW_file <- X[,c(1,2,6)]

#getting the node weights and matching to ew
("Matching gene-level p-values...")
shared.genes <- intersect(c(EW_file[,1],EW_file[,2]),node_weights[,1])
shared.nodes.weight <- node_weights[match(shared.genes, node_weights[,1]), ]
qnorm(shared.nodes.weight[,2]/2, lower.tail = F) -> node.z
shared.nodes.weight$z <- node.z

NW_file <- shared.nodes.weight[,c(1,3)]
write.table(NW_file, "NodeWeights.txt", row.names = F, col.names = F, quote = F, sep="\t")

shared.edge.weights <- X[which(X[,1] %in% shared.genes & X[,2] %in% shared.genes), ]
EW_file <- shared.edge.weights[,c(1,2,6)]
write.table(EW_file, "EdgeWeights.txt", row.names = F, col.names = F, quote = F, sep="\t")


message("Edge weight and node weight files have been created :)")