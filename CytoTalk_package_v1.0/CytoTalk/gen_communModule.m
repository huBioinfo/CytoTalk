%--Summarize all PCSFs.
clear
gen_summaryPCSF

%--Select the most robust PCSF from the background PCSFs.
clear
gen_robustPCSF

%--Generate files ready for import into Cytoscape.
clear
load pcsf_finalResult4 beta_final omega_final
gen_CytoFileFun(beta_final,omega_final);


exit  %-this is to exit MATLAB in order to run subsequent programs.


