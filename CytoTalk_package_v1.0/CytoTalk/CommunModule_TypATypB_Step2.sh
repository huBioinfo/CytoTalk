#!/bin/bash

#$ -o Step2.log

#$ -cwd


#4) Construct the integrated gene network.
date +%D-%H:%M:%S
echo "(4/6) Constructing the integrated gene network...(around 20 min)"
Rscript comp_NonSelfTalkScore_TypA.R
Rscript comp_NonSelfTalkScore_TypB.R
matlab -nodisplay -nodesktop -nosplash -r construct_integratedNetwork  #one should specify the absolute path of the executable Matlab program.
date +%D-%H:%M:%S

#5) Generate multiple PCSFs based on the integrated gene network. ~40min on 2 million edges.
echo "(5/6) Generating multiple PCSFs...(around 40 min)"
bash JobRun_Parallel.sh
date +%D-%H:%M:%S

#6) Generate the communication gene module.
echo "(6/6) Generating the final communication gene module...(around 5 min)"
matlab -nodisplay -nodesktop -nosplash -r gen_communModule  #one should specify the absolute path of the executable Matlab program.
echo "ALL DONE! Cytotalk has generated output files in the folder: /Cytotalk/IllustratePCSF/"
date +%D-%H:%M:%S


