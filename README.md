# CytoTalk<br />
### Version 1.2 (April 15th, 2020)<br />

### Overview<br />
Cell-cell communication in a tissue microenvironment is mediated by signal transduction pathways. Single-cell technology has opened the door for studying signal transduction at much higher resolution in a complex tissue. Currently, there is a lack of analytical methods to infer signaling pathways based on single-cell omics data. Here we introduce a computational method, **CytoTalk**, for _de novo_ construction of **cell type-specific signaling networks** using single-cell transcriptomics data. Using an integrated intracellular and intercellular gene network as the input, CytoTalk identifies candidate pathways using prize-collecting Steiner forest (PCSF) algorithm. CytoTalk is implemented using MATLAB (version >= R2018a), R (version >= 3.5.0) and Python (version >= 3.7.0).

<br />

<div align=center><img src="https://github.com/huBioinfo/CytoTalk/blob/master/CytoTalk_schematic.png" width="60%" height="60%" /></div>



### Required R and Python packages and system environment variables<br />
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
>>Tip: set the number of logical cores available in Line 13 of the R scripts “comp_MIcoexp_TypA_WinPara.R” and “comp_MIcoexp_TypB_WinPara.R”. The default is 14 logical cores. 

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








