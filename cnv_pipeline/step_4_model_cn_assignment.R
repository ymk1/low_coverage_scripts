
# Multi-segmentation Extract Segments and Estimates
mseg_gr_list<-list()

for(y in colnames(mcols(mseg_comp))[3:ncol(mcols(mseg_comp))]){
    
   # Extract and Reformat Sample
        s<-mseg_comp[,y]
        t.m<-split(s, values(s)[,1])
  
        t<-as.data.frame(reduce(t.m, 
                              min.gapwidth=100000L))
        colnames(t)[2]<-"sample"
        t<-t[,c(3:7,2)]
        t$sample<-2^as.numeric(t$sample)
  

    # Density-based Clustering
        d_clust<-dbscan::dbscan(as.matrix(as.numeric(t$sample)), 
                              eps = .1, minPts=2)
        CN2_group<-names(sort(-table(d_clust$cluster)))[1]
  
  
    # Assign Largest Group as CN2 and Obtain Median Est
        CN2<-median(as.numeric(t$sample)[which(d_clust$cluster==CN2_group)], 
                    na.rm = TRUE)
       # Note: Above does not work when assumption majority of genome is CN2
        if(CN2<0.825){
            CN2<-1
         }
    
    # Assign Second Largest Group (Lower than CN2) to be CN1 with Various Conditions:
      if(max(d_clust$cluster)>1){
              CN_2nd_peak<-names(sort(-table(d_clust$cluster)))[2]
              CN_2nd_peak<-abs(abs(median(as.numeric(temp$sample)[which(d_clust$cluster==CN_2nd_peak)], 
                                  na.rm = TRUE))-CN2)
      
        if(CN_2nd_peak<0.2 & max(d_clust$cluster)>2){
              CN_2nd_peak<-names(sort(-table(d_clust$cluster)))[3]
              CN_2nd_peak<-abs(median(as.numeric(temp$sample)[which(d_clust$cluster==CN_2nd_peak)], 
                                  na.rm = TRUE))
        if(CN_2nd_peak>0.7){
              CN_2nd_peak<-0.5
        }
      } 
        if(CN_2nd_peak<0.2 & max(d_clust$cluster)==2){
              CN_2nd_peak<-0.5
          }
      }
  
     if(max(d_clust$cluster)==1){
            CN_2nd_peak<-0.5
      }
  
     # Defining States and Distances
        d = sapply(c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6), 
                   function(x){abs(CN2)+abs(CN_2nd_peak)*x})
        alpha = rep(1, length(d))
        alpha = rbind(alpha/sum(alpha))

        j<-as.numeric(t$sample)
        post <- as.data.frame(sapply(d, function(rc) 
                      exp(abs(rc)*log(rc)-rc-lgamma(abs(j)+1))))
      for(z in 1:ncol(ll)){
              post[,z]<-post[,z]*alpha[z]}
        colnames(post)<-c("CN0", "CN0.5", "CN1", 
                          "CN1.5", "CN2", "CN2.5", 
                          "CN3", "CN4", "CN5",
                          "CN6", "CN7", "CN8")
  
        t$CN_call<-colnames(post)[apply(post, 1, which.max)]
        t$sample_name<-y
  
        t$state<-ifelse(t$CN_call %in% class[1:4], "loss",
                        ifelse(t$CN_call %in% "CN2", "neutral",
                            "gain"))
                                
  mseg_gr_list[[y]]<-t
  
}
