#Preparing RNA seq data

RNAseq <- read.table("C:/users/amanuel1/Dev/ZhaoLab/MS/GSE111972_norm_data.txt", header=T, as.is = T)
row.names(RNAseq) <- RNAseq[,1]
RNAseq <- RNAseq[,-c(1)]
RNAseq <- data.matrix(RNAseq)
RNA_t <- (quantile_normalisation(log2(RNAseq+1)))

GE.control <- RNA_t[,c(seq(1,16))]
GE.case <- RNA_t[,c(seq(17,31))]



GE.control.WM <- GE.control[,c(2,4,6,8,10,11,12,13,14,15,16)]
GE.case.WM <- GE.case.WM <- GE.case[,c(2,4,6,8,10,11,12,13,14,15)]


