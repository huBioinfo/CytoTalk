%-Combine analysis results in each folder together to record Edge Background.
clear
load pcsf_finalResult4 qualityTable
PCSF_edgeIndex_Tab=[];
for i=1:size(qualityTable,1)
    beta=qualityTable(i,1);
    omega=qualityTable(i,2);
    
    %---------progress bar-------------%
%     fprintf('beta %.6f, omega %.5f.\n',beta,omega); %only show two digits behind the dot.
    %----------------------------------%
    ResultFileName=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega),'/pcsf_result4');
        
    %--combine analysis results.    
    load(ResultFileName);
        
    PCSF_edgeIndex_Tab=[PCSF_edgeIndex_Tab;PCSF_edgeIndex];
    
    clearvars -except beta omega PCSF_edgeIndex_Tab qualityTable
end
save pcsf_finalResult4 PCSF_edgeIndex_Tab -append


%-Count the repeated times of each edge in each PCSF module.
clear
load pcsf_finalResult4 PCSF_edgeIndex_Tab qualityTable
PCSF_edgeCount_fullList=cell(size(qualityTable,1),1);
PCSF_edgeCount_medianList=zeros(size(PCSF_edgeCount_fullList));
for i=1:size(qualityTable,1)
    beta=qualityTable(i,1);
    omega=qualityTable(i,2);
    
    %---------progress bar-------------%
    fprintf('beta %.6f, omega %.5f.\n',beta,omega); %only show two digits behind the dot.
    %----------------------------------%
    ResultFileName=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega),'/pcsf_result4');
    
    %--combine analysis results.    
    load(ResultFileName);
        
    if ~isempty(PCSF_edgeIndex)
        edgeCount=zeros(length(PCSF_edgeIndex),1);
        for k=1:length(PCSF_edgeIndex)
            edgeCount(k)=length(find(PCSF_edgeIndex_Tab==PCSF_edgeIndex(k))); %-at least one edge.
        end
        
        PCSF_edgeCount_fullList{i}=edgeCount; 
        PCSF_edgeCount_medianList(i)=median(edgeCount);
    end
        
    clearvars -except beta omega PCSF_edgeIndex_Tab qualityTable PCSF_edgeCount_fullList PCSF_edgeCount_medianList
     
end
save pcsf_finalResult4 PCSF_edgeCount_fullList PCSF_edgeCount_medianList -append

clear
load pcsf_finalResult4 PCSF_edgeCount_fullList qualityTable
PCSF_edgeKSPval_fullList=ones(size(PCSF_edgeCount_fullList,1),3);
PCSF_edgeRankSumPval_fullList=ones(size(PCSF_edgeCount_fullList,1),3);
for i=1:size(qualityTable,1)
    beta=qualityTable(i,1);
    omega=qualityTable(i,2);
    
    %---------progress bar-------------%
    fprintf('beta %.6f, omega %.5f.\n',beta,omega); %only show two digits behind the dot.
    %----------------------------------%
    PCSF_edgeKSPval_fullList(i,1)=beta;
    PCSF_edgeKSPval_fullList(i,2)=omega;
    PCSF_edgeRankSumPval_fullList(i,1)=beta;
    PCSF_edgeRankSumPval_fullList(i,2)=omega;
    
    if ~isempty(PCSF_edgeCount_fullList{i})
        
        modifiedList=PCSF_edgeCount_fullList;
        %-remove self.
        modifiedList(i,:)=[];
        
        %-remove empty entries (no edges in the PCSF, only isolated root nodes after artificial edges removed).
        emptyInd=find(cellfun(@isempty,modifiedList));
        modifiedList(emptyInd,:)=[];
        
        %-compare to combined other parameter settings.
        other_edgeCount=[];
        for kk=1:length(modifiedList)
            other_edgeCount=[other_edgeCount;modifiedList{kk}];
        end
            
        %-using KS test.
        [h,PCSF_edgeKSPval_fullList(i,3)]=kstest2(PCSF_edgeCount_fullList{i},other_edgeCount,'Tail','smaller');

        %-using rank-sum test.
        PCSF_edgeRankSumPval_fullList(i,3)=ranksum(PCSF_edgeCount_fullList{i},other_edgeCount,'tail','right');
             
    end %-end of if.

end
save pcsf_finalResult4 PCSF_edgeKSPval_fullList PCSF_edgeRankSumPval_fullList -append

%-Select the PCSF with the minimum ks-pvalue.
clear
load pcsf_finalResult4 PCSF_edgeKSPval_fullList
[aa,indind]=sort(PCSF_edgeKSPval_fullList(:,3),'ascend');
beta_final=PCSF_edgeKSPval_fullList(indind(1),1);
omega_final=PCSF_edgeKSPval_fullList(indind(1),2);
save pcsf_finalResult4 beta_final omega_final -append


