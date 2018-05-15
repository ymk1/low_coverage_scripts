
library(BioStrings)
library(plyr)
library(aCGH)

#-------------------------------------------------------#
# Skipped Winsorization as PCA works, moving to Segmentation:

# Step 1: Prepare Data Frame for Segmentation
    tumors.pcf<-cbind(tumor[,1:2], tumors_pca)
    tumors.pcf[,1]<-gsub("Chr", "", tumors.pcf[,1])
    colnames(tumors.pcf)[1:2]<-c("chrom", "pos")
    chrom.arm<-rep('p', nrow(tumors.pcf))


# Step 2: Checking the Appropriate Gamma
    copynumber::plotGamma(tumors.pcf[,c(1,2,10)], cex = 3)

    
# Step 3: Split by Group and run Multipcf on each:
    #haplos <- read.csv("~/Documents/Tasmanian_Devil/sample database/2018-1-30_v21_Haplotype_Sample_Info_consolidated_EPM.csv", stringsAsFactors=FALSE)
    unique_haplos<-plyr::count(haplos$Groups)
 # Run Multiseg based on Haplogroups:
    multiseg_list<-list()
    multiseg_estimates<-list()

  for(x in which(unique_haplos[,2]>1)){
  multiseg<-copynumber::multipcf(tumors.pcf[,c(1,2,which(gsub("X", "", colnames(tumors.pcf)) %in% 
                                                     haplos$Sample.TCG_ID[which(haplos$Groups %in% 
                                                                                  unique_haplos[x,1])]))], 
                                 arms=chrom.arm, 
                                 gamma = 50,
                                 return.est = TRUE)
  multiseg_list[[unique_haplos[x,1]]]<-as.data.frame(multiseg$segments)
  multiseg_estimates[[unique_haplos[x,1]]]<-as.data.frame(multiseg$estimates)
  
  }
  multiseg_list<-multiseg_list[unlist(lapply(multiseg_list, length) != 0)]
  multiseg_estimates<-multiseg_estimates[unlist(lapply(multiseg_estimates, length) != 0)]

  # Compile Multisegmentation:
  multiseg_compiled<-Reduce(function(x, y) merge(x, y, all = TRUE), 
                            multiseg_estimates)
  multiseg_compiled<-multiseg_compiled[order(multiseg_compiled$pos),]
  multiseg_compiled<-multiseg_compiled[order(multiseg_compiled$chrom),]


for(y in unique(colnames(multiseg_compiled)[3:ncol(multiseg_compiled)])){
  message(y)
  if(!is.na(match(y, colnames(multiseg_compiled)))){
  multiseg_compiled[,y]<-mergeLevels(tumors.pcf[, y], multiseg_compiled[,y], verbose = FALSE)$vecMerged
}
}


multiseg_compiled$end<-multiseg_compiled$pos+99999
colnames(multiseg_compiled)[1:2]<-c("seqnames", "start")

multiseg_compiled<-makeGRangesFromDataFrame(multiseg_compiled, keep.extra.columns = TRUE)
