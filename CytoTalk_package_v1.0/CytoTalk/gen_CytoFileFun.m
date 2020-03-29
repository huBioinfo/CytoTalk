function [] = gen_CytoFileFun(beta_final,omega_final)
%GEN_CYTOFILEFUN Summary of this function goes here
%   Detailed explanation goes here

ResultFileName=strcat('./bt',num2str(beta_final,'%.6f'),'/bt',num2str(beta_final,'%.6f'),'_omg',num2str(omega_final),'/pcsf_result4');
load(ResultFileName);

load Exp_cleaned_2 TypAExp_rmRep_dmNew
load Exp_cleaned_4 TypBExp_rmRep_dmNew
load IntegratedNet_TypATypB_ID

%---------Extract gene expression of resulting PCSF.
PCSF_geneExp=zeros(length(PCSF_nodeSym),1);
PCSF_geneCellType=zeros(length(PCSF_nodeSym),1); %-"1" as TypA; "2" as TypB.
for i=1:length(PCSF_nodeSym)
    %---------progress bar-------------%
%     fprintf('PCSF_nodeSym %d.\n',i);
    %----------------------------------%
    if strcmp(PCSF_nodeSym{i}(end-5:end),'__TypA') %-if true, it's TypA.
        PCSF_geneCellType(i)=1; %-"1" as TypA.
        vec=strcmp(PCSF_nodeSym{i}(1:end-6),TypAExp_rmRep_dmNew.RowNames); 
        ind=find(vec); %-must have a single value.
        PCSF_geneExp(i)=mean(double(TypAExp_rmRep_dmNew(ind,:)));
    elseif strcmp(PCSF_nodeSym{i}(end-5:end),'__TypB') %-if true, it's TypB.
        PCSF_geneCellType(i)=2; %-"2" as TypB.
        vec=strcmp(PCSF_nodeSym{i}(1:end-6),TypBExp_rmRep_dmNew.RowNames); 
        ind=find(vec); %-must have a single value.
        PCSF_geneExp(i)=mean(double(TypBExp_rmRep_dmNew(ind,:)));
    end
end

%---------Extract gene prize of resulting PCSF.
PCSF_genePrize=zeros(length(PCSF_nodeSym),1);
for i=1:length(PCSF_nodeSym)
    %---------progress bar-------------%
%     fprintf('PCSF_nodeSym %d.\n',i); 
    %----------------------------------%
    vec=strcmp(PCSF_nodeSym(i),integratedNet_GeneSym);
    ind=find(vec);
    PCSF_genePrize(i)=integratedNet_GenePrize_initial(ind);
end

%---------Extract gene Real Name of resulting PCSF.
PCSF_geneRealName=cell(length(PCSF_nodeSym),1);
for i=1:length(PCSF_nodeSym)
    %---------progress bar-------------%
%     fprintf('PCSF_nodeSym %d.\n',i); 
    %----------------------------------%
    PCSF_geneRealName{i}=PCSF_nodeSym{i}(1:end-6);
end
save pcsf_resultCyto PCSF_nodeSym PCSF_geneExp PCSF_geneCellType PCSF_edgeSym PCSF_edgeIndex PCSF_edgeCost PCSF_edgeCellType isolatedNode PCSF_genePrize PCSF_geneRealName


%-------------------Generate files for Cytoscape--------------------------%
%-create a folder for Cytoscape files.
status_1=mkdir('./IllustratePCSF/');

fid1=fopen('./IllustratePCSF/PCSF_edgeSym.sif','w');
formatSpec_1='%s\t%s\t%s\n';
for i=1:size(PCSF_edgeSym,1)
    
    fprintf(fid1,formatSpec_1,PCSF_edgeSym{i,1},'pp',PCSF_edgeSym{i,2});
    
end
fclose(fid1);

fid2=fopen('./IllustratePCSF/PCSF_edgeCost.txt','w');
formatSpec_2='%s %s %s\t%f\n';
fprintf(fid2,'%s\t%s\n','EdgeSym','EdgeCost');
for i=1:size(PCSF_edgeSym,1)

    fprintf(fid2,formatSpec_2,PCSF_edgeSym{i,1},'(pp)',PCSF_edgeSym{i,2},PCSF_edgeCost(i));
    
end
fclose(fid2);

fid5=fopen('./IllustratePCSF/PCSF_edgeCellType.txt','w');
formatSpec_5='%s %s %s\t%d\n';
fprintf(fid5,'%s\t%s\n','EdgeSym','EdgeCellType');
for i=1:size(PCSF_edgeSym,1)

    fprintf(fid5,formatSpec_5,PCSF_edgeSym{i,1},'(pp)',PCSF_edgeSym{i,2},PCSF_edgeCellType(i));
    
end
fclose(fid5);

fid3=fopen('./IllustratePCSF/PCSF_geneExp.txt','w');
formatSpec_3='%s\t%f\n';
fprintf(fid3,'%s\t%s\n','GeneSym','GeneExp');
for i=1:length(PCSF_nodeSym)

    fprintf(fid3,formatSpec_3,PCSF_nodeSym{i},PCSF_geneExp(i));
    
end
fclose(fid3);

fid4=fopen('./IllustratePCSF/PCSF_geneCellType.txt','w');
formatSpec_4='%s\t%d\n';
fprintf(fid4,'%s\t%s\n','GeneSym','GeneCellType');
for i=1:length(PCSF_nodeSym)

    fprintf(fid4,formatSpec_4,PCSF_nodeSym{i},PCSF_geneCellType(i));
    
end
fclose(fid4);

fid6=fopen('./IllustratePCSF/PCSF_genePrize.txt','w');
formatSpec_6='%s\t%f\n';
fprintf(fid6,'%s\t%s\n','GeneSym','GenePrize');
for i=1:length(PCSF_nodeSym)

    fprintf(fid6,formatSpec_6,PCSF_nodeSym{i},PCSF_genePrize(i));
    
end
fclose(fid6);

fid7=fopen('./IllustratePCSF/PCSF_geneRealName.txt','w');
formatSpec_7='%s\t%s\n';
fprintf(fid7,'%s\t%s\n','GeneSym','GeneRealName');
for i=1:length(PCSF_nodeSym)
    
    fprintf(fid7,formatSpec_7,PCSF_nodeSym{i},PCSF_geneRealName{i});
    
end
fclose(fid7);


%-----------Generate a txt file for crosstalk edges along with normalized crosstalk score--------%
%-prepocess the format of crosstalk edge symbol.
for i=1:size(crosstalk_edgeSym,1)
    if strcmp(crosstalk_edgeSym{i,1}(end-5:end),'__TypB') %-if true, it's TypB, then switch the two symbols.
        c=crosstalk_edgeSym{i,1};
        crosstalk_edgeSym{i,1}=crosstalk_edgeSym{i,2};
        crosstalk_edgeSym{i,2}=c;
    end
end
%-sort the formatted crosstalk edges based on normalized crostalk score.
[crosstalk_edgeScoSorted,indind]=sort(crosstalk_edgeSco,'descend');
crosstalk_edgeSymSorted=crosstalk_edgeSym(indind,:);

fid8=fopen('./IllustratePCSF/PCSF_CrosstalkEdge.txt','w');
formatSpec_8='%s\t%s\t%f\n';
fprintf(fid8,'%s\t%s\t%s\n','Gene_in_CellTypeA','Gene_in_CellTypeB','NormalizedCrosstalkScore');
for i=1:size(crosstalk_edgeSymSorted,1)

    fprintf(fid8,formatSpec_8,crosstalk_edgeSymSorted{i,1}(1:end-6),crosstalk_edgeSymSorted{i,2}(1:end-6),crosstalk_edgeScoSorted(i));
    
end
fclose(fid8);
        

end


