%--Process two intracellular networks based on matrices generated in the Step1.
clear
gen_2intracellularNet

%-----Compute node prize using PEM score---------------------%
clear
comp_NodePrize_CellType

clear
comp_NodePrize_RWR  %~14min

clear
comp_NodePrize_Combine
%-------------------Gene Node Prize Done!!!-----------------------%

%------Compute crosstalk score based on non-self-talk score---------%
clear
comp_Crosstalk_specific
%-----------------------Crosstalk Score Done!!!---------------------%

%-------------Below for network analysis--------------------------%
clear
comp_topNet_TypA

clear
comp_topNet_TypATypB

clear
comp_CrosstalkNet

clear
comp_ArtiNet  %-artiNode is connected to all genes of the integrated networks.

clear
comp_IntegratedNet

clear
write_IntegratedNet

clear
gen_diffOmgBt_Parallel_v2  %6min,10G/615run-Parallelly generate different folders for different parameter (omega and beta) settings for parallel PCSF running.

exit  %-this is to exit MATLAB in order to run subsequent programs.


