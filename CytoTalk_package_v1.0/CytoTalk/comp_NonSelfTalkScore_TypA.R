##--compute NonSelfTalk score of TypA part for all known ligand-receptor pairs.
library(entropy)

NonSelfTalkSym_Tumo <- vector()  #-initiated as empty.
NonSelfTalkSco_Tumo <- vector()

species <- read.table("Species.txt", header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

if (species$V1 == "Human") {
  
  LRpair <- read.table("LigandReceptor_Human.txt", header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
}

if (species$V1 == "Mouse") {
  
  LRpair <- read.table("LigandReceptor_Mouse.txt", header = FALSE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
}


#load Tumo expression matrix
TumoExp0 <- read.csv("scRNAseq_CellTypeA.csv") #now row names as first column.
TumoExp <- TumoExp0[, -1] #remove the first column.
rownames(TumoExp) <- TumoExp0[, 1]  #set first column as row names.

totalNum <- 1:nrow(LRpair)
for (i in totalNum) {
  #---------progress bar-------------#
  #cat("LRpair:", i, "\n");
  #----------------------------------#
  
  ind_1 <- match(LRpair[i, 1], row.names(TumoExp))
  
  if (!is.na(ind_1)) {
    ind_2 <- match(LRpair[i, 2], row.names(TumoExp))
    
    if (!is.na(ind_2)) {
      
          Tumo_exp_g1 <- TumoExp[ind_1,]  #-checked.
          Tumo_exp_g2 <- TumoExp[ind_2,] 
          
          #-number of bins is set to sqrt(number of samples)--WGCNA mannual.
          #-compute mutual information in Tumo.
          
          #-Here need to add noise to all-zero-genes.
          Tumo_exp_g1_value <- as.numeric(Tumo_exp_g1)
          Tumo_exp_g2_value <- as.numeric(Tumo_exp_g2)
          
          if (all(Tumo_exp_g1_value==0)) {
            randCell_1 <- sample(1:length(Tumo_exp_g1_value), 1)
            Tumo_exp_g1_value[randCell_1] <- 1E-20
          }
          
          if (all(Tumo_exp_g2_value==0)) {
            randCell_2 <- sample(1:length(Tumo_exp_g2_value), 1)
            Tumo_exp_g2_value[randCell_2] <- 1E-20
          }
          
          Tumo_table <- discretize2d(Tumo_exp_g1_value, Tumo_exp_g2_value, numBins1 = sqrt(ncol(Tumo_exp_g1)), numBins2 = sqrt(ncol(Tumo_exp_g2)))
          Tumo_H1 <- entropy(rowSums(Tumo_table), method = "MM")  #compute marginal entropies.
          Tumo_H2 <- entropy(colSums(Tumo_table), method = "MM")
          Tumo_H12 <- entropy(Tumo_table, method = "MM")
          Tumo_MI <- Tumo_H1+Tumo_H2-Tumo_H12
          
          normFactor <- min(c(Tumo_H1,Tumo_H2))
          Tumo_MI_dist <- -log10(Tumo_MI/normFactor)
 
          NonSelfTalkSym_Tumo <- rbind(NonSelfTalkSym_Tumo, LRpair[i,])
          NonSelfTalkSco_Tumo <- rbind(NonSelfTalkSco_Tumo, Tumo_MI_dist) 
          
    }
  }
  
}  #-end of For loop.

write.table(NonSelfTalkSym_Tumo, "./NonSelfTalkSym_TypA.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)   #-checked.
write.table(NonSelfTalkSco_Tumo, "./NonSelfTalkSco_TypA.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)   #-checked.


