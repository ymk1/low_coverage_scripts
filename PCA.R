# Internal Function
log2R_samples<-function(x){
  b<-log2(x/median(x, na.rm = TRUE))
  b[which(is.infinite(b))]<-0
  return(b)
}

# Note: this is not the GATK version:
# DataFrames Format: Chrom, Start, End, Host name 1, Host name 2, ...)

PCA_normalization<-function(host_df, tumor_df){

# Step 1: LogR the Raw Read Counts
host_log<-apply(host[,4:ncol(host)], 2, log2R_samples)
tumor_log<-apply(tumor[,c(4:ncol(tumor))], 2, log2R_samples)


# Step 2: Calculate PCA
res<-prcomp(t(host_log), center = TRUE, scale = FALSE)

# Selected the First 4 PC:
test1<-predict(res, newdata=t(tumor_log))[,1:4] %*% t(res$rotation[,1:4])

# Created your noise signal for each tumor and centered it on your PON
trunc <- scale(test1, center = -1 * res$center, scale=FALSE)

# Remove the noise signal from your tumor
tumors_pca<-sapply(c(1:ncol(tumor_log)), function(x){
  df<-tumor_log[,x]-trunc[x,]
})
colnames(tumors_pca)<-colnames(tumor_log)

return(tumors_pca)}
