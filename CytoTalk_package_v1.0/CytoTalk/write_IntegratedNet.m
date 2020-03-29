%-------------------Write edgeID and edge costs and node prizes (later) of integrated network to txt file for identifying PCSF using Python.
clear
load IntegratedNet_TypATypB_ID
fid5=fopen('IntegratedNet_edge.txt','w');
formatSpec_1='%d\t%d\n'; %this is tab-delimited format for importing to Python.
for pp=1:size(integratedNet_EdgeID,1)
    %---------progress bar-------------%
    %fprintf('Edge %d.\n',pp);
    %----------------------------------%
    fprintf(fid5,formatSpec_1,integratedNet_EdgeID(pp,1),integratedNet_EdgeID(pp,2));
end
fclose(fid5);

%-replace -0 with 0 before writing it to txt file. %-use the non-fixed edge costs.
zeroInd=find(integratedNet_EdgeCost_common==-0); 
integratedNet_EdgeCost_common(zeroInd)=0; %-replacing -0 with 0.
fid6=fopen('IntegratedNet_edgeCost_common.txt','w');
formatSpec_2='%f\n'; 
for qq=1:length(integratedNet_EdgeCost_common)
    %---------progress bar-------------%
    %fprintf('Edge %d.\n',qq);
    %----------------------------------%
    fprintf(fid6,formatSpec_2,integratedNet_EdgeCost_common(qq));
end
fclose(fid6);

%-write RootNode.txt
fid8=fopen('RootNode.txt','w');
formatSpec_3='%d\n'; 
fprintf(fid8,formatSpec_3,length(integratedNet_GeneSym)-1);
fclose(fid8);


