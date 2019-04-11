#Preparing RNA seq data

RNAseq <- read.table("C:/users/amanuel1/Dev/ZhaoLab/MS/GSE111972_norm_data.txt", header=T, as.is = T)
row.names(RNAseq) <- RNAseq[,1]
RNAseq <- RNAseq[,-c(1)]
RNAseq <- data.matrix(RNAseq)


GE.control <- RNAseq[,c(seq(1,16))]
GE.case <- RNAseq[,c(seq(17,31))]

GE.control.WM <- GE.control[,c(2,4,6,8,10,11,12,13,14,15,16)]
GE.case.WM <- GE.case.WM <- GE.case[,c(2,4,6,8,10,11,12,13,14,15)]


