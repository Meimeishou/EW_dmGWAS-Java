
#Calculating edge weights for EW_dmGWAS input, matching PPI edges to node weights.
#04/10/2019

GE.control <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CONTROL")
GE.case <- read.table("PATH_TO_GENE_EXPRESSION_DATA_FOR_CASE")
PPI = read.table("PATH_TO_PPI")


#Obtaining PCC
obtainPCC <- function(GE.matrix){
  PCCvalues <- cor(t(GE.matrix))
  return(PCCvalues)
}

controlPCC <- obtainPCC(GE.control)
casePCC <- obtainPCC(GE.case)

#matching to PPI
matchPCCtoPPI <- function(PCCvalues, PPI){
  
  #match to PPI
  return(PPI_PCC_rmNA_rmSI)
}


control_PPI_PCC <- matchPCCtoPPI(controlPCC, PPI)
case_PPI_PCC <- matchPCCtoPPI(casePCC, PPI)


