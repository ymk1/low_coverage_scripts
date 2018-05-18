# Orthogonal Approach to Possion:
# Credit to Kevin


# Internal Functions:
# ASCAT function
get_r <- function(rho, psi, n) {
  log2((2 * (1-rho) + rho * n) / psi)
}

# Least Squares Sum Calculation
compare <- function(observed, ideal) {
  sq <- (observed - ideal)^2
  sum(sq)
}

# Comparison
objective <- function(x) {
  rho <- x[1]
  psi <- x[2]
  n = 2
  ideal <- get_r(rho, psi, n)
  compare(observed, ideal)
}

# CN State
cn_call <- function(x, ct, p){
  cn <- (2*exp(x)-p*(1-ct))/ct
  round(cn)
 }

# Libraries:
library(nloptr)

#------------------------#

# Currently Running Across Entire Genome, but may be more optimal to run on Density Based Defined CN2 cluster?
# HAVE NOT OPTIMIZED FOR ACROSS ALL SAMPLES
         
# Per Sample Basis, need tumor cont fraction and raw read counts.
ct = haplos$Tumor.Purity[which(haplos$Sample.TCG_ID %in% y)]
observed = tumors.pcf[-grep("X", tumors.pcf[,1]),y]


# Optimize Parameters
opt <- nloptr( x0=c(ct,2,2), # ploidy default: 2; CN default: 2
                eval_f=objective,
                lb = c(0,0,0),
                ub = c(1,5,4),
                opts = list("algorithm"="NLOPT_LN_COBYLA",
                            "xtol_rel"=1.0e-8))

pars<-opt$solution
names(pars)<-c("ct", "ploidy", "cn")


# Incorporate Parameters to Estimate CN States from Segment Output
s.sseg<-sseg.gr[which(sseg.gr$sampleID %in% y),]
s.sseg$CN<-cn_call(sseg.gr$mean, pars[1], pars[2])
         
s.mseg<-mseg.gr[which(mseg.gr$sampleID %in% y),]
s.mseg$CN<-cn_call(mseg.gr$mean, pars[1], pars[2])

  
  














