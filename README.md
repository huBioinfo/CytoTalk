CytoTalk 
================
**R Version 4.0.0 (October 7th, 2021)**

**Please refer to https://github.com/tanlabcode/CytoTalk for the latest release, which is a completely re-written R version of the program. Approximately half of the run time as been shaved off, the program is now cross-compatible with Windows and \*NIX systems, the file space usage is down to roughly a tenth of what it was, and graphical outputs have been made easier to import or now produce portable SVG files with embedded hyperlinks.**

R Version 3.1.0 (June 8th, 2021)

R Version 3.0.3 (May 31st, 2021)

R Version 3.0.2 (May 19th, 2021)

R Version 3.0.1 (May 12th, 2021)

MATLAB Version 2.1 (see Releases page)

MATLAB Version 2.0 (February 22nd, 2021) 


## Overview

Signal transduction is the primary mechanism for cell-cell communication. scRNA-Seq technology holds great promise for studying cell-cell communication at much higher resolution. Signaling pathways are highly dynamic and cross-talk among them is prevalent. Due to these two features, simply examining expression levels of ligand and receptor genes cannot reliably capture the overall activities of signaling pathways and interactions among them. 

We have developed the **CytoTalk algorithm for *de novo* construction of a signaling network (union of multiple signaling pathways emanating from the ligand-receptor pairs) between two cell types using single-cell transcriptomics data.** The algorithm first constructs an integrated network consisting of intracellular and intercellular functional gene interactions. It then identifies the signaling network by solving a prize-collecting Steiner forest (PCSF) problem based on appropriately defined node prize (i.e. cell-specific gene activity) and edge cost (i.e. probability of functional interaction between two genes). The objective of the PCSF problem is to **find an optimal subnetwork in the integrated network that includes genes with high levels of cell-type-specific expression and close connection to highly active ligand-receptor pairs**. 
<br />
<br />

⚠ **Important usage tips for the old version v3.1.0 (Please go to https://github.com/tanlabcode/CytoTalk for the latest CytoTalk_v4.0.0)**

* Download "CytoTalk_package_v3.1.0.zip" from the Master branch or the Releases page (https://github.com/huBioinfo/CytoTalk/releases/tag/v3.1.0) and refer to the user manual inside the package.

* Gene expression matrix files (named as "scRNAseq_xxx.csv") for ALL cell types in the tissue or tumor microenvironment under study are required to REPLACE 7 example "scRNAseq_BCells/EndothelialCells/TCells/xxx.csv" in the Input/ folder, in order to compute cell-type-specificity of gene expression built in the CytoTalk algorithm. Please note that:

    (1) The row names of all these gene expression matrices should be exactly the same;

    (2) The gene expression values should be ln-transformed normalized scRNA-seq data (e.g. Seurat-preprocessed data);

    (3) The Input/ folder should ONLY contain preprocessed "scRNAseq_xxx.csv" mentioned above and several existing .txt, .R and .py files, no other .csv files.
<br />
<br />

⚠ **Update log** 

**2021-10-07: The latest release “CytoTalk_v4.0.0” is a completely re-written R version of the program. Approximately half of the run time as been shaved off, the program is now cross-compatible with Windows and \*NIX systems, the file space usage is down to roughly a tenth of what it was, and graphical outputs have been made easier to import or now produce portable SVG files with embedded hyperlinks.**

2021-06-08: The release "CytoTalk_v3.1.0" is a major updated R version on the basis of v3.0.3. We have added a function to generate Cytoscape files for visualization of each ligand-receptor-associated pathway extracted from the predicted signaling network between the two given cell types. For each predicted ligand-receptor pair, its associated pathway is defined as the user-specified order of the neighborhood of the ligand and receptor in the two cell types.

2021-05-31: The release "CytoTalk_v3.0.3" is a revised R version on the basis of v3.0.2. A bug has been fixed in this version to avoid errors occurred in some special cases. We also provided a new example "RunCytoTalk_Example_StepByStep.R" to run the CytoTalk algorithm in a step-by-step fashion. Please download "CytoTalk_package_v3.0.3.zip" from the Releases page (https://github.com/huBioinfo/CytoTalk/releases/tag/v3.0.3) and refer to the user manual inside the package.

2021-05-19: The release "CytoTalk_v3.0.2" is a revised R version on the basis of v3.0.1. A bug has been fixed in this version to avoid running errors in some extreme cases. Final prediction results will be the same as v3.0.1. Please download the package from the Releases page (https://github.com/huBioinfo/CytoTalk/releases/tag/v3.0.2) and refer to the user manual inside the package.

2021-05-12: The release "CytoTalk_v3.0.1" is an R version, which is more easily and friendly to use!! Please download the package from the Releases page (https://github.com/huBioinfo/CytoTalk/releases/tag/v3.0.1) and refer to the user manual inside the package.




<div align=center><img src="https://github.com/huBioinfo/CytoTalk/blob/master/CytoTalk_schematic.png" width="60%" height="60%" /></div>
<br />


     
## Cite CytoTalk  

* Hu Y, Peng T, Gao L, Tan K. CytoTalk: *De novo* construction of signal transduction networks using single-cell transcriptomic data. ***Science Advances***, 2021, 7(16): eabf1356.

    https://advances.sciencemag.org/content/7/16/eabf1356

* Hu Y, Peng T, Gao L, Tan K. CytoTalk: *De novo* construction of signal transduction networks using single-cell RNA-Seq data. *bioRxiv*, 2020.
 
    https://www.biorxiv.org/content/10.1101/2020.03.29.014464v1

## Reference  

* Shannon P, et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Research*, 2003, 13: 2498-2504.

## Contact  

Kai Tan, tank1@chop.edu<br />
Yuxuan Hu, huyuxuan@xidian.edu.cn<br />


