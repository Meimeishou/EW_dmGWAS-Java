Edge weight
# Normalization
# Handle dupliated probes
# Calculate the PCC
# Match with PPI network
# Define new statistic
## Based on the equations in the EW_dmGWAS paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514922/
> df$F_r_case <- log((1+ df$case_cor)/(1-df$case_cor))/2
> df$F_r_control <- log((1+ df$control_cor)/(1-df$control_cor))/2
> df$X <- (df$F_r_case - df$F_r_control)/ sqrt(1/4+1/3)
> df$e <- qnorm(1-2*(1-pnorm(abs(df$X))))
> dim(df)
[1] 33829     8
> df <- na.omit(df)
> dim(df)
[1] 31855     8
> scale(df$e) -> df$e_scale
> head(df)

> write.table(df[,c(1,2,9)],"edge_weight.txt",sep="\t",row.names=FALSE)

