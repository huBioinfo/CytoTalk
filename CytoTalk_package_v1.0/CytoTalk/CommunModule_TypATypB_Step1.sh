#!/bin/bash

#$ -o Step1.log

#$ -cwd


#Identify the communication module between Cell Type A and Cell Type B.
#Input: "scRNAseq_CellTypeA/B/C/D/E/.../X.csv", and "Species.txt".
#Output: An output folder, /Cytotalk/IllustratePCSF/, contains a network file "PCSF_edgeSym.sif" and six attribute files that are ready for import into Cytoscape.

#1) Pre-process the scRNA-Seq data.
echo "(1/6) Pre-processing the scRNA-Seq data...(very fast)"
matlab -nodisplay -nodesktop -nosplash -r gen_intracellularNetMat  #one should specify the absolute path of the executable Matlab program.

#2) Generate mutual information matrix using parallel environment.~6-12h depending on the number of single cells and available logical cores.
date +%D-%H:%M:%S
echo "(2/6) Computing mutual information for cell type A...(around 4 hours)"
Rscript comp_MIcoexp_TypA_WinPara.R
date +%D-%H:%M:%S
echo "(2/6) Computing mutual information for cell type B...(around 4 hours)"
Rscript comp_MIcoexp_TypB_WinPara.R

#3) Generate indirect edge-filtered network matrix using parallel environment.~30min
date +%D-%H:%M:%S
echo "(3/6) Removing indirect edges for cell type A...(around 15 min)"
Rscript comp_GeneNet_TypA_LinuxPara.R
date +%D-%H:%M:%S
echo "(3/6) Removing indirect edges for cell type B...(around 15 min)"
Rscript comp_GeneNet_TypB_LinuxPara.R
echo "Now it's ready for Step2 to identify the communication gene module."
date +%D-%H:%M:%S


