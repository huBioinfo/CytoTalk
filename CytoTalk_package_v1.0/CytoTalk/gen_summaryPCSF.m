%--Copy analysis files.
clear
%-First beta range.
for beta=5:5:300  %-modify 1/3 for cancer study.
    %-Second omega range.
    for omega=0.1:0.1:1.5
        %---------progress bar-------------%
%         fprintf('beta %.6f, omega %.5f.\n',beta,omega); %only show two digits behind the dot.
        %----------------------------------%
        BtOmgFolderName=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega),'/');
        
        %--copy PCSF result analysis files.    
        str8_source='deter_W4.m';
        status_8=copyfile(str8_source,BtOmgFolderName);
        
%         str9_source='comp_FisherMethodPval.R';
%         status_9=copyfile(str9_source,BtOmgFolderName);
     end %-end of omega.
     
end  %-end of beta.

%-Run analysis files in each folder.
clear
%-First beta range.
for beta=5:5:300  %-modify 2/3 for cancer study.
    %-Second omega range.
    for omega=0.1:0.1:1.5
        %---------progress bar-------------%
        fprintf('beta %.6f, omega %.5f.\n',beta,omega); %only show two digits behind the dot.
        %----------------------------------%
        BtOmgFolderName=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega),'/');
        
        %--run analysis files.    
        oldName=cd(BtOmgFolderName);
        
        deter_W4;
        cd ../../
     end %-end of omega.
     
end  %-end of beta.


%-Quality check for the generated PCSFs.
clear
qualityTable=[];
for beta=5:5:300   %-modify 3/3 for cancer study.The upper bound is not fixed,must above 2000 edges.
    
    for omega=0.1:0.1:1.5   %-This range is fixed due to very consistent result using different datasets.
        %---------progress bar-------------%
%         fprintf('beta %.6f, omega %.5f.\n',beta,omega); %only show two digits behind the dot.
        %----------------------------------%
        ResultFileName=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega),'/pcsf_result4');
        
        %--combine analysis results.    
        load(ResultFileName);    
        qualityList=gen_qualityCheck_v2(beta,omega,PCSF_nodeSym,PCSF_edgeSym,crosstalk_edgeSym,nTree,isolatedNode);
        
        qualityTable=[qualityTable;qualityList];
        
        clearvars -except beta omega qualityTable

     end %-end of omega.
     
end  %-end of beta.

%-clean data_1---remove beta values that make the number of edges in the PCSF
%above 2000 edges, which is our empirical cutoff.
rmRowInd_beta=find(qualityTable(:,15)>=2000); %Column15 is the number of all edges in the PCSF.
qualityTable(rmRowInd_beta,:)=[];

%-clean data_2---remove omega values==0.1/0.2/0.3/0.4 due to very small
%number of PCSF edges.
rmRowInd_omega=find(qualityTable(:,2)<0.5);  %Column2 is the Omega value.
qualityTable(rmRowInd_omega,:)=[];
save pcsf_finalResult4 qualityTable


