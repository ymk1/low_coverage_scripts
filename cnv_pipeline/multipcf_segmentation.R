
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
    mseg_list<-list()
    mseg_est<-list()

  for(x in which(unique_haplos[,2]>1)){
        # Run Multipcf on A Haplogroup      
        multiseg<-copynumber::multipcf(tumors.pcf[,c(1,2,which(colnames(tumors.pcf) %in% haplos$Sample.TCG_ID[which(haplos$Groups %in% 
                                                                                  unique_haplos[x,1])]))], 
                                        arms=chrom.arm, 
                                        gamma = 50,
                                        return.est = TRUE)
      
        # Add Data to List
        #mseg_list[[unique_haplos[x,1]]]<-as.data.frame(mseg$segments)
        mseg_est[[unique_haplos[x,1]]]<-as.data.frame(mseg$estimates)
  
  }
        #mseg_list<-mseg_list[unlist(lapply(mseg_list, length) != 0)]
        mseg_est<-mseg_est[unlist(lapply(mseg_est, length) != 0)]


# Step 4: Compile Multisegmentation across All Haplogroups and Reduce & Sort Matrix
        mseg_comp<-Reduce(function(x, y) merge(x, y, all = TRUE), 
                            mseg_est)
        mseg_comp<-mseg_comp[order(mseg_comp$pos),]
        mseg_comp<-mseg_comp[order(mseg_comp$chrom),]


# Step 5: Perform Ansari/Bradley Tests (pval<.05) for Preliminary Merging for Each Sample
for(y in unique(colnames(mseg_comp)[3:ncol(mseg_comp)])){
    
    message(y)
    
        if(!is.na(match(y, colnames(mseg_comp)))){
            mseg_comp[,y]<-mergeLevels(tumors.pcf[, y], 
                                               mseg_comp[,y], 
                                               verbose = FALSE)$vecMerged
            }
    }


# Step 6: Reformat Data Frames for CN calling assignments
    mseg_comp$end<-mseg_comp$pos+99999
    colnames(mseg_comp)[1:2]<-c("seqnames", "start")

    mseg_comp<-makeGRangesFromDataFrame(mseg_comp, 
                                        keep.extra.columns = TRUE)
