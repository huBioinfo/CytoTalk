%----Import TypA mutual info and indirect edge-filtered network----------------%
clear
load Exp_cleaned_2 TypAExp_rmRep_dmNew
nGene=size(TypAExp_rmRep_dmNew,1);
fid4=fopen('GeneNetwork_TypA.txt');
fmt_1 = ['%*s' repmat(' %f',1,nGene)];
MiMatrix_TypA=textscan(fid4,fmt_1,'delimiter',' ','Headerlines',1);
MiMatrix_TypA=cell2mat(MiMatrix_TypA);
fclose(fid4);
save MI_TypA MiMatrix_TypA -v7.3

clear
load Exp_cleaned_2 TypAExp_rmRep_dmNew
RowNames_TypA=TypAExp_rmRep_dmNew.RowNames;
save MI_TypA RowNames_TypA -append

%-make it to list format.
clear
load MI_TypA 
isolatedNode_fromRow=[];
isolatedNode_fromCol=[];
for i=1:size(MiMatrix_TypA,1)
    if isempty(find(MiMatrix_TypA(i,:)))
        isolatedNode_fromRow=[isolatedNode_fromRow;i];
    end
end
for j=1:size(MiMatrix_TypA,2)
    if isempty(find(MiMatrix_TypA(:,j)))
        isolatedNode_fromCol=[isolatedNode_fromCol;j];
    end
end
if isequal(isolatedNode_fromRow,isolatedNode_fromCol)
    MiMatrix_TypA(isolatedNode_fromRow,:)=[];
    MiMatrix_TypA(:,isolatedNode_fromCol)=[];  %-make sure the matrix is symmetric.
    RowNames_TypA(isolatedNode_fromRow,:)=[];
end

if length(unique(RowNames_TypA))==length(RowNames_TypA)  %-make sure RowNames is unique.
    correctFlag=1;
else
    correctFlag=0;
end

MiList_value_TypA=zeros(nchoosek(size(RowNames_TypA,1),2),1);
MiList_genePair_TypA=zeros(nchoosek(size(RowNames_TypA,1),2),2);
k=1;
for i=1:(size(RowNames_TypA,1)-1)
    %---------progress bar-------------%
    %fprintf('RowNames %d.\n',i);
    %----------------------------------%
    for j=(i+1):size(RowNames_TypA,1)
        if MiMatrix_TypA(i,j)~=0
            MiList_genePair_TypA(k,1)=i;
            MiList_genePair_TypA(k,2)=j;
            MiList_value_TypA(k)=MiMatrix_TypA(i,j);
            k=k+1;
        end
    end
end 
MiList_value_TypA(k:length(MiList_value_TypA),:)=[];
MiList_genePair_TypA(k:size(MiList_genePair_TypA,1),:)=[];
MiList_genePair_TypA=sort(MiList_genePair_TypA,2); %-sort each row for convenient comparison in the future.
save MI_TypA MiList_value_TypA MiList_genePair_TypA isolatedNode_fromRow isolatedNode_fromCol correctFlag -append
%-----------------------TypA done!

%----Import TypB mutual info and indirect edge-filtered network---------------%
clear
load Exp_cleaned_4 TypBExp_rmRep_dmNew
nGene=size(TypBExp_rmRep_dmNew,1);
fid4=fopen('GeneNetwork_TypB.txt');
fmt_1 = ['%*s' repmat(' %f',1,nGene)];
MiMatrix_TypB=textscan(fid4,fmt_1,'delimiter',' ','Headerlines',1);
MiMatrix_TypB=cell2mat(MiMatrix_TypB);
fclose(fid4);
save MI_TypB MiMatrix_TypB -v7.3

clear
load Exp_cleaned_4 TypBExp_rmRep_dmNew
RowNames_TypB=TypBExp_rmRep_dmNew.RowNames;
save MI_TypB RowNames_TypB -append

%-make it to list format.
clear
load MI_TypB  
isolatedNode_fromRow=[];
isolatedNode_fromCol=[];
for i=1:size(MiMatrix_TypB,1)
    if isempty(find(MiMatrix_TypB(i,:)))
        isolatedNode_fromRow=[isolatedNode_fromRow;i];
    end
end
for j=1:size(MiMatrix_TypB,2)
    if isempty(find(MiMatrix_TypB(:,j)))
        isolatedNode_fromCol=[isolatedNode_fromCol;j];
    end
end
if isequal(isolatedNode_fromRow,isolatedNode_fromCol)
    MiMatrix_TypB(isolatedNode_fromRow,:)=[];
    MiMatrix_TypB(:,isolatedNode_fromCol)=[];  %-make sure the matrix is symmetric.
    RowNames_TypB(isolatedNode_fromRow,:)=[];
end

if length(unique(RowNames_TypB))==length(RowNames_TypB)  %-make sure RowNames is unique.
    correctFlag=1;
else
    correctFlag=0;
end

MiList_value_TypB=zeros(nchoosek(size(RowNames_TypB,1),2),1);
MiList_genePair_TypB=zeros(nchoosek(size(RowNames_TypB,1),2),2);
k=1;
for i=1:(size(RowNames_TypB,1)-1)
    %---------progress bar-------------%
    %fprintf('RowNames %d.\n',i);
    %----------------------------------%
    for j=(i+1):size(RowNames_TypB,1)
        if MiMatrix_TypB(i,j)~=0
            MiList_genePair_TypB(k,1)=i;
            MiList_genePair_TypB(k,2)=j;
            MiList_value_TypB(k)=MiMatrix_TypB(i,j);
            k=k+1;
        end
    end
end   %-checked!
MiList_value_TypB(k:length(MiList_value_TypB),:)=[];
MiList_genePair_TypB(k:size(MiList_genePair_TypB,1),:)=[];
MiList_genePair_TypB=sort(MiList_genePair_TypB,2); %-sort each row for convenient comparison in the future.
save MI_TypB MiList_value_TypB MiList_genePair_TypB isolatedNode_fromRow isolatedNode_fromCol correctFlag -append
%-------------------------TypB done!


