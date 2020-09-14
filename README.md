# CytoTalk<br />
### Version 1.2 (April 15th, 2020)<br />
<br />

## Overview<br />
Cell-cell communication in a tissue microenvironment is mediated by signal transduction pathways. Single-cell technology has opened the door for studying signal transduction at much higher resolution in a complex tissue. Currently, there is a lack of analytical methods to infer signaling pathways based on single-cell omics data. Here we introduce a computational method, **CytoTalk**, for _de novo_ construction of **cell type-specific signaling networks** using single-cell transcriptomics data. Using an integrated intracellular and intercellular gene network as the input, CytoTalk identifies candidate pathways using prize-collecting Steiner forest (PCSF) algorithm. CytoTalk is implemented using MATLAB (version >= R2018a), R (version >= 3.5.0) and Python (version >= 3.7.0).

<br />

<div align=center><img src="https://github.com/huBioinfo/CytoTalk/blob/master/CytoTalk_schematic.png" width="60%" height="60%" /></div>
<br />

## Prepare input files<br />
1. CytoTalk requires a comma-delimited “.csv” file containing scRNA-Seq data for each cell type under study. Each file contains the log2-transformed normalized scRNA-Seq data for a cell type with rows as genes (GENE SYMBOL) and columns as cells.<br />

>The files should be named as:<br />
>**“scRNAseq_CellTypeA.csv”<br />
>“scRNAseq_CellTypeB.csv”<br />
>“scRNAseq_CellTypeC.csv”<br />
>“scRNAseq_CellTypeD.csv”<br />
>“scRNAseq_CellTypeE.csv”<br />
>…<br />**
>>Tips:
>>>(1) Csv files for all cell types should be copied into the /CytoTalk/ folder since all of them are needed for computing cell-type-specificity of gene expression in the CytoTalk algorithm.

>>>(2) The /CytoTalk/ folder can only be used **once** for predicting the signaling network between cell type A and cell type B. Please make sure the gene expression data of your interested cell type pair are stored in “scRNAseq_CellTypeA.csv” and “scRNAseq_CellTypeB.csv”. Examples are in the /ExampleInput/ folder.

2. CytoTalk also requires a **“Species.txt"** file indicating the species from which the scRNA-Seq data are generated. Currently, “Human” and “Mouse” are supported. An example is in the folder /ExampleInput/ folder. This file should be also copied into the /CytoTalk/ folder.
<br />

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
Yuxuan Hu, yuxuan_hu_xd@163.com<br />


