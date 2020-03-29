clear
%-find the file names of all cell types.
allExpFile=dir('./scRNAseq_CellType*.csv');
allExpVector=cell(1,length(allExpFile));
for i=1:length(allExpFile)
    fileName=allExpFile(i).name;
    allExpVector{i}=bioma.data.DataMatrix('File',fileName,'Delimiter',',');
end
save Exp_allCSV allExpFile allExpVector -v7.3

clear
load Exp_allCSV allExpVector allExpFile
allExpVector_NoLog=cell(size(allExpVector));
%-need to convert back to count-based data. Currently, Exp=log2(count+1).
for i=1:length(allExpVector)
    ExpNoLog_rmRep_dm=2.^allExpVector{i}; %-checked.
    allExpVector_NoLog{i}=ExpNoLog_rmRep_dm-1;   %-checked.
end
save Exp_allCSV_NoLog allExpVector_NoLog allExpFile -v7.3

clear
load Exp_allCSV_NoLog allExpVector_NoLog allExpFile
same_rowname=allExpVector_NoLog{1}.RowNames;
save GeneCellTypeSpecific same_rowname allExpFile

clear
load Exp_allCSV_NoLog allExpVector_NoLog
Exp_tpmMean=cell(size(allExpVector_NoLog));
Exp_tpmSum=cell(size(allExpVector_NoLog));
for i=1:length(allExpVector_NoLog)
    Exp_tpm=double(allExpVector_NoLog{i});
    Exp_tpmMean{i}=mean(Exp_tpm,2);
    Exp_tpmSum{i}=sum(Exp_tpmMean{i}); 
end
save GeneCellTypeSpecific Exp_tpmMean Exp_tpmSum -append

%--Compute gene summary and dataset summary used for PEM score.
clear
load GeneCellTypeSpecific Exp_tpmMean Exp_tpmSum
geneMatrix=[];
datasetSum=0;
for i=1:length(Exp_tpmMean)
    geneMatrix=[geneMatrix,Exp_tpmMean{i}];
    datasetSum=datasetSum+Exp_tpmSum{i};
end
geneSum=sum(geneMatrix,2);
save GeneCellTypeSpecific geneSum datasetSum -append
    
clear
load GeneCellTypeSpecific same_rowname geneSum datasetSum Exp_tpmMean Exp_tpmSum
typeSpecific=cell(size(Exp_tpmMean));
for kk=1:length(Exp_tpmMean)
    typeSpecific{kk}=zeros(length(same_rowname),1);
    for i=1:length(same_rowname)
        e=geneSum(i)*(Exp_tpmSum{kk}/datasetSum);
        typeSpecific{kk}(i)=log10(Exp_tpmMean{kk}(i)/e);
    end
end
save GeneCellTypeSpecific typeSpecific -append


