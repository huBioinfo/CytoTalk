library(infotheo)
library(doParallel)

#load Macr expression matrix
MacrExp0 <- read.csv("TypBExp_rmRepGf_dm.csv") #now row names as first column.
MacrExp <- MacrExp0[, -1] #remove the first column.
rownames(MacrExp) <- MacrExp0[, 1]  #set first column as row names.

MacrExp <- t(MacrExp)
MacrExp_discrete <- discretize(MacrExp)

#-setup parallel environment.
registerDoParallel(cores=6) #-number of physical cores. Can also be the number of logical cores (current=6).
#getDoParWorkers()

#calculate mutual information matrix.
#Sys.time()
MImatrix <- foreach (i = 1:ncol(MacrExp_discrete), .packages = c('infotheo', 'doParallel'), .combine = 'rbind') %dopar% { 
  
  #---------progress bar-------------#
  #cat("Gene:", i, "\n");
  #----------------------------------#
  
  foreach (j = 1:ncol(MacrExp_discrete), .combine = 'c') %do% {
    
    mutinformation(MacrExp_discrete[, i], MacrExp_discrete[, j], method = "mm") 
    
  }
}
#Sys.time()

write.table(MImatrix, "./MutualInfo_TypB_Para.txt", quote = FALSE)
#Sys.time()


