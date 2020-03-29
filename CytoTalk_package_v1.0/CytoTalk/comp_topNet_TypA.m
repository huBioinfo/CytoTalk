%--------------------Generate edge costs of TypA for PCSF.
clear
load MI_TypA
mu=mean(MiList_value_TypA);
sigma=std(MiList_value_TypA);
MiList_value_TypA=(MiList_value_TypA-mu)/sigma;
save MI_topNet_TypA MiList_genePair_TypA MiList_value_TypA

clear
load MI_topNet_TypA MiList_genePair_TypA
load MI_TypA RowNames_TypA
GeneInEdge=unique([MiList_genePair_TypA(:,1);MiList_genePair_TypA(:,2)]);
if length(GeneInEdge)==length(RowNames_TypA)  %-should hold, if not, matrix to list part is wrong.
    MiList_geneSym_TypA=cell(length(RowNames_TypA),1);
    for pp=1:length(RowNames_TypA)
        %---------progress bar-------------%
        %fprintf('Gene %d.\n',pp);
        %----------------------------------%
        MiList_geneSym_TypA{pp}=strcat(RowNames_TypA{pp},'__TypA');
    end
end
save MI_topNet_TypA MiList_geneSym_TypA -append

        
%--------------------Generate node prize of TypA for PCSF.
clear
load MI_topNet_TypA MiList_geneSym_TypA
load GeneNodePrize GenePrize_TypA RowNames_TypA
MiList_genePrize_TypA=zeros(length(MiList_geneSym_TypA),1);
for i=1:length(MiList_geneSym_TypA)
    %---------progress bar-------------%
    %fprintf('geneSym %d.\n',i);
    %----------------------------------%
    vec=strcmp(MiList_geneSym_TypA{i}(1:end-6),RowNames_TypA); %-extract original gene symbol without "__TypA".
    ind=find(vec); %-must have a single value.
    MiList_genePrize_TypA(i)=GenePrize_TypA(ind);
end
save MI_topNet_TypA MiList_genePrize_TypA -append


