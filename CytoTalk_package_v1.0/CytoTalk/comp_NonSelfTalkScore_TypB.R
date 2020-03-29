##--compute NonSelfTalk score of TypB part for all known ligand-receptor pairs.
library(entropy)

NonSelfTalkSym_Macr <- vector()  #-initiated as empty.
NonSelfTalkSco_Macr <- vector()

species <- read.table("Species.txt", header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

if (species$V1 == "Human") {
  
  LRpair <- read.table("LigandReceptor_Human.txt", header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
}

if (species$V1 == "Mouse") {
  
  LRpair <- read.table("LigandReceptor_Mouse.txt", header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
}


#load Macr expression matrix
MacrExp0 <- read.csv("scRNAseq_CellTypeB.csv") #now row names as first column.
MacrExp <- MacrExp0[, -1] #remove the first column.
rownames(MacrExp) <- MacrExp0[, 1]  #set first column as row names.

totalNum <- 1:nrow(LRpair)
for (i in totalNum) {
  #---------progress bar-------------#
  #cat("LRpair:", i, "\n");
  #----------------------------------#
  
  ind_1 <- match(LRpair[i, 1], row.names(MacrExp))
  
  if (!is.na(ind_1)) {
    ind_2 <- match(LRpair[i, 2], row.names(MacrExp))
    
    if (!is.na(ind_2)) {
      
          Macr_exp_g1 <- MacrExp[ind_1,]  #-checked.
          Macr_exp_g2 <- MacrExp[ind_2,] 
          
          #-number of bins is set to sqrt(number of samples)--WGCNA mannual.
          #-compute mutual information in Macr.
          
          #-Here need to add noise to all-zero-genes.
          Macr_exp_g1_value <- as.numeric(Macr_exp_g1)
          Macr_exp_g2_value <- as.numeric(Macr_exp_g2)
          
          if (all(Macr_exp_g1_value==0)) {
            randCell_1 <- sample(1:length(Macr_exp_g1_value), 1)
            Macr_exp_g1_value[randCell_1] <- 1E-20
          }
          
          if (all(Macr_exp_g2_value==0)) {
            randCell_2 <- sample(1:length(Macr_exp_g2_value), 1)
            Macr_exp_g2_value[randCell_2] <- 1E-20
          }
          
          Macr_table <- discretize2d(Macr_exp_g1_value, Macr_exp_g2_value, numBins1 = sqrt(ncol(Macr_exp_g1)), numBins2 = sqrt(ncol(Macr_exp_g2)))
          Macr_H1 <- entropy(rowSums(Macr_table), method = "MM")  #compute marginal entropies.
          Macr_H2 <- entropy(colSums(Macr_table), method = "MM")
          Macr_H12 <- entropy(Macr_table, method = "MM")
          Macr_MI <- Macr_H1+Macr_H2-Macr_H12
          
          normFactor <- min(c(Macr_H1,Macr_H2))
          Macr_MI_dist <- -log10(Macr_MI/normFactor)
 
          NonSelfTalkSym_Macr <- rbind(NonSelfTalkSym_Macr, LRpair[i,])
          NonSelfTalkSco_Macr <- rbind(NonSelfTalkSco_Macr, Macr_MI_dist) 
          
    }
  }
  
}  #-end of For loop.

write.table(NonSelfTalkSym_Macr, "./NonSelfTalkSym_TypB.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)   #-checked.
write.table(NonSelfTalkSco_Macr, "./NonSelfTalkSco_TypB.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)   #-checked.


