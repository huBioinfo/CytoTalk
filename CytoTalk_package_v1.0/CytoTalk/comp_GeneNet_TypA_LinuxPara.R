#load Tumo mutual information matrix.
MIadj_Tumo <- read.table("MutualInfo_TypA_Para.txt", header = TRUE, sep = " ") #-checked.

#set diagonal values to zeros to avoid self-cycles.
totalGene <- 1:nrow(MIadj_Tumo)
for (i in totalGene) {
  #---------progress bar-------------#
  #cat("Gene:", i, "\n"); #Very fast.
  #----------------------------------#
  
  MIadj_Tumo[i,i] <- 0
  
}

#remove indirect edges.(should set environment variable: export OMP_NUM_THREADS=n)
library(parmigene)
#Sys.time()
MIadj_Tumo <- as.matrix(MIadj_Tumo)
GN_Tumo <- aracne.m(MIadj_Tumo, 0.15) 
#Sys.time()

write.table(GN_Tumo, "./GeneNetwork_TypA.txt", quote = FALSE) 
#Sys.time()


