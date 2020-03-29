#load Macr mutual information matrix.
MIadj_Macr <- read.table("MutualInfo_TypB_Para.txt", header = TRUE, sep = " ") #-checked.

#set diagonal values to zeros to avoid self-cycles.
totalGene <- 1:nrow(MIadj_Macr)
for (i in totalGene) {
  #---------progress bar-------------#
  #cat("Gene:", i, "\n"); #Very fast.
  #----------------------------------#
  
  MIadj_Macr[i,i] <- 0
  
}

#remove indirect edges.(should set environment variable: export OMP_NUM_THREADS=n)
library(parmigene)
#Sys.time()
MIadj_Macr <- as.matrix(MIadj_Macr)
GN_Macr <- aracne.m(MIadj_Macr, 0.15) 
#Sys.time()

write.table(GN_Macr, "./GeneNetwork_TypB.txt", quote = FALSE) 
#Sys.time()


