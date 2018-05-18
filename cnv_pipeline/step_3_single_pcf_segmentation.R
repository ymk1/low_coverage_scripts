
library(BioStrings)
library(GRanges)
library(copynumber)
library(aCGH)

#-------------------------------------------------#
# Call Single Sample Segmentation across the entire dataset:

# Step 1: Run Single-pcf
     sseg<-copynumber::pcf(tumors.pcf, kmin = 10,
                            arms=chrom.arm, 
                            gamma = 70, 
                            return.est = TRUE)


# Step 2: Extract Segments and Estimates from the Single-segmentation 
     sseg.gr<-as.data.frame(sseg$segments)[,c(2,4,5,6,1,7)]
     sseg.gr_est<-as.data.frame(sseg$estimates)
    
    
# Step 3: Perform Ansari/Bradley Tests (pval<.05) for Preliminary Merging for Each Sample
 for(y in unique(sseg.gr$sampleID)){
      
          message(y)
  
          sseg.gr_est[,y]<-mergeLevels(tumors.pcf[, y], 
                                          sseg.gr_est[,y], 
                                          verbose = FALSE)$vecMerged
      }


# Step 4: Convert the Data to GRanges Format
    sseg.gr_cest<-with(sseg.gr_est[,c(1,2,which(colnames(sseg.gr_est)==y))], 
                                    GRanges(seqnames = chrom,
                                    IRanges(start = pos, 
                                            end = pos+100000-1)))
                                                                                      
 for(y in unique(sseg.gr$sampleID)){
    
           mcols(sseg.gr_cest)$mean<-sseg.gr_est[,y]
           colnames(mcols(sseg.gr_cest))[ncol(mcols(sseg.gr_cest))]<-y
    }

    sseg.gr_cest<-as.data.frame(sseg.gr_cest)
    colnames(sseg.gr_cest)<-gsub("X", "", colnames(single_seg_cn_estimates))
    
    
    


