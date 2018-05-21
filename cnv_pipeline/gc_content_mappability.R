# Correct using the correct.GC function (optimized further using the HMM Copy to include mappability)
# GC Correction Function
correct.GC <- function(rc.f, # read count data (each row is an interval)
                       gc.df=GC.content,
                       map.df=map.gc.content,
                       samplesize = 50000){ # Calculated GC Fraction for each window bin
  
  # a. Create a GC Factor Class for the LOESS Regression
        gc.class = cut(gc.df$GCcontent, # Cut the GC Fraction (0-1) into 0.02 intervals
                       breaks = seq(0, 1, 0.02), 
                       include.lowest = TRUE)

  
  # b. Assign each window read count to one of the GC Factor Classes
        samp.ii = unlist(tapply(1:length(gc.class), gc.class, 
                                function(e) e[sample(1:length(e), 
                                                     min(c(length(e), 500)))]))


  # c. Add the Contig and Position Information to the data frame of read counts   
        rc.df<-data.frame(gc.df[,1:3], rc.f)
        colnames(rc.df)<-c("chr", "start", "end", "rc") 
        rc.df$GCcontent = gc.df$GCcontent
        rc.df$MPcontent = map.df$Mapcontent

  
  # d. Calculate the LOESS Regression and residual      
    # create a regression based on the read count and GC content interval
        lo = stats::loess(rc ~ GCcontent, data = rc.df[samp.ii, ])

    # calculate the residual from the regression to the raw
        rc.df$cor.gc = mean(rc.df$rc, na.rm = TRUE) * rc.df$rc/stats::predict(lo, newdata = rc.df)

  # e. Mappability using LOESS Regression
    # i. Create a Map Factor Class for the LOESS Regression
        mp.class = cut(rc.df$MPcontent, # Cut the GC Fraction (0-1) into 0.02 intervals
                       breaks = seq(0, 1, 0.02), 
                       include.lowest = TRUE)
  
  # ii. Assign each window read count to one of the GC Factor Classes
        samp.ii = unlist(tapply(1:length(mp.class), mp.class, 
                                function(e) e[sample(1:length(e), 
                                                     min(c(length(e), 500)))]))

  # iii. Calculate the LOESS Regression and residual      
    # create a regression based on the read count and GC content interval
        lo = stats::loess(cor.gc ~ MPcontent, data = rc.df[samp.ii, ])
  
    # calculate the residual from the regression to the raw
        rc.df$cor.map = mean(rc.df$cor.gc, na.rm = TRUE) * rc.df$cor.gc/stats::predict(lo, newdata = rc.df)  
  
  # Mappability Correction
        #coutlier<-0.1
        #rc.df$valid<-TRUE
        #rc.df$valid[ rc.df$rc <= 0 | rc.df$GCcontent < 0]<- FALSE

        #crange<-quantile(rc.df$cor.gc[rc.df$valid], 
        #             prob = c(0, 1-coutlier), na.rm = TRUE)
        #set<-which(rc.df$cor.gc<crange[2])
        #select<-sample(set, min(length(set), samplesize))
        #final = approxfun(lowess(rc.df$MPcontent[select],
        #                     rc.df$cor.gc[select]))

        #rc.df$cor.map<- mean(rc.df$cor.gc, na.rm = TRUE) * rc.df$cor.gc/final(rc.df$MPcontent)      

  
  
        rc.df$cor.map = round(rc.df$cor.map, digits = 2)
        if (any(rc.df$cor.map < 0, na.rm = TRUE)) 
        rc.df$rc[rc.df$cor.map < 0] = 0
        rc.df$GCcontent = NULL

  return(rc.df$cor.map)
} 
