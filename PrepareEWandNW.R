
#Calculating edge weights for EW_dmGWAS input, matching PPI edges to node weights.
#04/10/2019

GE.control <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CONTROL", as.is = T)
GE.case <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CASE", as.is = T)
PPI <- read.table("PATH_TO_PPI")


#Obtaining PCC
obtainPCC <- function(GE.matrix){
  PCCvalues <- cor(t(GE.matrix))
  return(PCCvalues)
}

controlPCC <- obtainPCC(GE.control)
casePCC <- obtainPCC(GE.case)

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


control_PPI_PCC <- matchPCCtoPPI(controlPCC, PPI)
case_PPI_PCC <- matchPCCtoPPI(casePCC, PPI)


