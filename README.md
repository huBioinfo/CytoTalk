CytoTalk 
================
Version 2.0 (February 22nd, 2021) 


## Overview

* Signal transduction is the primary mechanism for cell-cell communication. scRNA-Seq technology holds great promise for studying cell-cell communication at much higher resolution. Signaling pathways are highly dynamic and cross-talk among them is prevalent. Due to these two features, simply examining expression levels of ligand and receptor genes cannot reliably capture the overall activities of signaling pathways and interactions among them. 

* We have developed the CytoTalk algorithm for de novo construction of a signaling network (union of multiple signaling pathways emanating from the ligand-receptor pairs) between two cell types using single-cell transcriptomics data. The algorithm first constructs an integrated network consisting of intracellular and intercellular functional gene interactions. It then identifies the signaling network by solving a prize-collecting Steiner forest (PCSF) problem based on appropriately defined node prize (i.e. cell-specific gene activity) and edge cost (i.e. probability of functional interaction between two genes). The objective of the PCSF problem is to find an optimal subnetwork in the integrated network that includes genes with high levels of cell-type-specific expression and close connection to highly active ligand-receptor pairs. CytoTalk is currently implemented using a combination of MATLAB (version  R2018a), R (version  3.5.0) and Python (version  3.7.0). 


<div align=center><img src="https://github.com/huBioinfo/CytoTalk/blob/master/CytoTalk_schematic.png" width="60%" height="60%" /></div>
<br />

## I. Input files 
(1) A comma-delimited “.csv” file containing scRNA-Seq data for each cell type under study. Each file contains the ln-transformed normalized scRNA-Seq data for a cell type with rows as genes (GENE SYMBOL) and columns as cells. Examples are in the /Input/ folder. The files should be named as:
“scRNAseq_Fibroblasts.csv”
“scRNAseq_Macrophages.csv”
“scRNAseq_EndothelialCells.csv”
“scRNAseq_CellTypeName.csv”
…

(2) A “TwoCellTypes.txt” file indicating the two cell types between which the signaling network is predicted. Please make sure that the cell type names should be consistent with scRNA-Seq data files above.

(3) A “Species.txt” file indicating the species from which the scRNA-Seq data are generated. Currently, “Human” and “Mouse” are supported.

(4) A “Cutoff_GeneFilter.txt” file indicating the cutoff for removing lowly-expressed genes in the processing of scRNA-Seq data. The default cutoff value is 0.1, which means that genes expressed in less than 10% of all cells of a given type are removed.

(5) A “BetaUpperLimit.txt” file indicating the upper limit of the test values of the algorithm parameter β, which is inversely proportional to the total number of genes in a given cell-type pair after removing lowly-expressed genes in the processing of scRNA-Seq data. Based on preliminary tests, the upper limit of β value is suggested to be 100 (default) if the total number of genes in a given cell-type pair is above 10,000. However, if the total number of genes is below 5000, it is necessary to increase the upper limit of β value to 500.

!!!Note that all example input files are in the /Input/ folder and should be customized and copied into the /CytoTalk/ folder before running. The /CytoTalk/ folder can only be used ONCE for a given cell-type pair. Please use a new /CytoTalk/ folder for analysis of other cell-type pairs.
## Install required R and Python packages and set system environment variables<br />
1. The following four R packages should be installed (R version ≥ 3.5.0 is recommended). 

>>(1) “entropy”: https://cran.r-project.org/web/packages/entropy/index.html
```R
install.packages("entropy") #R
```

>>(2) “infotheo”: https://cran.r-project.org/web/packages/infotheo/index.html
```R
install.packages("infotheo") #R
```

>>(3) “doParallel”: https://cran.r-project.org/web/packages/doParallel/index.html
```R
install.packages("doParallel") #R
```
>>Tip: Set the number of logical cores available in Line 13 of the R scripts “comp_MIcoexp_TypA_WinPara.R” and “comp_MIcoexp_TypB_WinPara.R”. The default is 14 logical cores. 

>>(4) “parmigene”: https://cran.rstudio.com/web/packages/parmigene/index.html
```R
install.packages("parmigene") #R
```
```Bash
export OMP_NUM_THREADS=n #Bash; n is the number of logical cores available.
```

2. The following three Python package should be installed (Python version ≥ 3.7.0 is recommended).

>>(1) “pcsf_fast”: https://github.com/fraenkel-lab/pcst_fast
```Bash
pip3 install pcst_fast #Bash
```
```Bash
export PYTHONPATH=$PYTHONPATH:/your installed pcsf_fast folder/
```

>>(2) “numpy”: https://numpy.org/install/
```Bash
pip3 install numpy #Bash
```

>>(3) “datetime”: https://pypi.org/project/DateTime/
```Bash
pip3 install DateTime #Bash
```

3. Set system environment variable to include the absolute path of the executable MATLAB program. An example in the macOS system is as following:
```Bash
export PATH=/Applications/MATLAB_R2018a.app/bin/:$PATH #Bash
```
<br />

## Run CytoTalk<br />
Copy the input file-added “/CytoTalk/” folder to your working directory and execute the following two steps:
```Bash
bash CommunModule_TypATypB_Step1.sh #Bash
```
>>Tip: This step may take up to 6 hours because computing mutual information for all gene pairs is time consuming.

```Bash
bash CommunModule_TypATypB_Step2.sh #Bash
```
<br />

## CytoTalk output
The output folder, “/CytoTalk/IllustratePCSF/”, contains a network file **“PCSF_edgeSym.sif”** and the following six attribute files that are ready for import into Cytoscape for visualization and further analysis of the predicted signaling network between cell type A and cell type B.
<br />

>Two edge attribute files:

>**(1) “PCSF_edgeCellType.txt”<br />
>(2) “PCSF_edgeCost.txt”<br />**

>Four node attribute files:

>**(1) “PCSF_geneCellType.txt”<br />
>(2) “PCSF_geneExp.txt”<br />
>(3) “PCSF_genePrize.txt”<br />
>(4) “PCSF_geneRealName.txt”<br />**
<br />

## Cite CytoTalk:
Hu, Y., Peng, T., Gao, L., & Tan, K. (2020). CytoTalk: _De novo_ construction of signal transduction networks using single-cell RNA-Seq data. _bioRxiv_.<br />
<br />

## Contact:
Kai Tan, tank1@email.chop.edu<br />


