%--------------------Generate edge costs of TypB for PCSF.
clear
load MI_TypB
mu=mean(MiList_value_TypB);
sigma=std(MiList_value_TypB);
MiList_value_TypB=(MiList_value_TypB-mu)/sigma;
save MI_topNet_TypB MiList_genePair_TypB MiList_value_TypB

clear
load MI_topNet_TypB MiList_genePair_TypB MiList_value_TypB
load MI_TypB RowNames_TypB
load MI_topNet_TypA MiList_genePair_TypA MiList_value_TypA MiList_geneSym_TypA
GeneInEdge=unique([MiList_genePair_TypB(:,1);MiList_genePair_TypB(:,2)]);
if length(GeneInEdge)==length(RowNames_TypB)  %-should hold, if not, matrix to list part is wrong.
    MiList_geneSym_TypB=cell(length(RowNames_TypB),1);
    for pp=1:length(RowNames_TypB)
        %---------progress bar-------------%
        %fprintf('Gene %d.\n',pp);
        %----------------------------------%
        MiList_geneSym_TypB{pp}=strcat(RowNames_TypB{pp},'__TypB');
    end
end

MiList_geneSym_TypATypB=[MiList_geneSym_TypA;MiList_geneSym_TypB];
MiList_genePair_TypB_new=MiList_genePair_TypB+length(MiList_geneSym_TypA);
MiList_genePair_TypATypB=[MiList_genePair_TypA;MiList_genePair_TypB_new];
MiList_value_TypATypB=[MiList_value_TypA;MiList_value_TypB];
save MI_topNet_TypB MiList_geneSym_TypB MiList_genePair_TypB_new -append
save MI_topNet_TypATypB MiList_geneSym_TypATypB MiList_genePair_TypATypB MiList_value_TypATypB

        
%--------------------Generate node prize of TypB for PCSF.
clear
load MI_topNet_TypB MiList_geneSym_TypB
load GeneNodePrize GenePrize_TypB RowNames_TypB
MiList_genePrize_TypB=zeros(length(MiList_geneSym_TypB),1);
for i=1:length(MiList_geneSym_TypB)
    %---------progress bar-------------%
    %fprintf('geneSym %d.\n',i);
    %----------------------------------%
    vec=strcmp(MiList_geneSym_TypB{i}(1:end-6),RowNames_TypB); %-extract original gene symbol without "__TypB".
    ind=find(vec); %-must have a single value.
    MiList_genePrize_TypB(i)=GenePrize_TypB(ind);
end
save MI_topNet_TypB MiList_genePrize_TypB -append

clear
load MI_topNet_TypA MiList_genePrize_TypA
load MI_topNet_TypB MiList_genePrize_TypB
MiList_genePrize_TypATypB=[MiList_genePrize_TypA;MiList_genePrize_TypB];
save MI_topNet_TypATypB MiList_genePrize_TypATypB -append


