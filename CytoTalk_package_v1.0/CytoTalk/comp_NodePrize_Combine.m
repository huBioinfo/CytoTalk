clear
load GeneCellTypeSpecific same_rowname typeSpecific allExpFile
load MI_TypA RowNames_TypA
load GeneRelevanceTypA Relevance_TypA
load MI_TypB RowNames_TypB
load GeneRelevanceTypB Relevance_TypB

for s=1:length(allExpFile)
    if strcmp('CellTypeA.csv',allExpFile(s).name(end-12:end))
       TypA_ind=s;
    end
    if strcmp('CellTypeB.csv',allExpFile(s).name(end-12:end))
       TypB_ind=s;
    end
end

TypeSpecific_TypA=zeros(length(RowNames_TypA),1);
for p=1:length(RowNames_TypA)
    ind=find(strcmp(RowNames_TypA(p),same_rowname));
    TypeSpecific_TypA(p)=typeSpecific{TypA_ind}(ind); %-must have only one value.
end

TypeSpecific_TypB=zeros(length(RowNames_TypB),1);
for q=1:length(RowNames_TypB)
    ind=find(strcmp(RowNames_TypB(q),same_rowname));
    TypeSpecific_TypB(q)=typeSpecific{TypB_ind}(ind); %-must have only one value.
end
save GeneNodePrize RowNames_TypA RowNames_TypB TypeSpecific_TypA TypeSpecific_TypB Relevance_TypA Relevance_TypB

clear
load GeneNodePrize
%-set under zero to 0 to make gene prize nonnegative.
underExpressInd_TypA=find(TypeSpecific_TypA<0);
TypeSpecific_TypA(underExpressInd_TypA)=0;
underExpressInd_TypB=find(TypeSpecific_TypB<0);
TypeSpecific_TypB(underExpressInd_TypB)=0;

GenePrize_TypA=Relevance_TypA.*TypeSpecific_TypA;
GenePrize_TypB=Relevance_TypB.*TypeSpecific_TypB;
save GeneNodePrize GenePrize_TypA GenePrize_TypB -append


