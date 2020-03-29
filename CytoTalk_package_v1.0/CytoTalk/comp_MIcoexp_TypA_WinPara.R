library(infotheo)
library(doParallel)

#load Tumo expression matrix
TumoExp0 <- read.csv("TypAExp_rmRepGf_dm.csv") #now row names as first column.
TumoExp <- TumoExp0[, -1] #remove the first column.
rownames(TumoExp) <- TumoExp0[, 1]  #set first column as row names.

TumoExp <- t(TumoExp)
TumoExp_discrete <- discretize(TumoExp)

#-setup parallel environment.
registerDoParallel(cores=6) #-number of physical cores. Can also be the number of logical cores (current=6).
#getDoParWorkers()

#calculate mutual information matrix.
#Sys.time()
MImatrix <- foreach (i = 1:ncol(TumoExp_discrete), .packages = c('infotheo', 'doParallel'), .combine = 'rbind') %dopar% { 
  
  #---------progress bar-------------#
  #cat("Gene:", i, "\n");
  #----------------------------------#
  
  foreach (j = 1:ncol(TumoExp_discrete), .combine = 'c') %do% {
    
    mutinformation(TumoExp_discrete[, i], TumoExp_discrete[, j], method = "mm") 
    
  }
}
#Sys.time()

write.table(MImatrix, "./MutualInfo_TypA_Para.txt", quote = FALSE)
#Sys.time()


