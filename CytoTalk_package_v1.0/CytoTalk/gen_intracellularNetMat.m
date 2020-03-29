%-1) Filter out lowly-expressed genes in cell type A.
%Data cleaning-1.
clear
TypAExp_rmRep_dm=bioma.data.DataMatrix('File', 'scRNAseq_CellTypeA.csv','Delimiter',',');
TypAExp_rmRep=double(TypAExp_rmRep_dm);
filterGene_TypA=[];
cutoff_TypA=floor(size(TypAExp_rmRep,2)*0.1);  %-must express in more than 10% (CellPhoneDB) of the cells in a specific cluster.
for p=1:size(TypAExp_rmRep,1)
    if length(find(TypAExp_rmRep(p,:)))<=cutoff_TypA
        filterGene_TypA=[filterGene_TypA;p];
    end
end
if ~isempty(filterGene_TypA)
    TypAExp_rmRep_dm(filterGene_TypA,:)=[];
end
save Exp_cleaned_1 TypAExp_rmRep_dm -v7.3

%-2) Only keep protein-coding genes.
%Data cleaning-2.
clear
load Exp_cleaned_1 TypAExp_rmRep_dm
fid0=fopen('Species.txt');
species=textscan(fid0,'%s','delimiter','\t');
fclose(fid0);
species=species{1};
if isequal(species{1},'Human')
    load GencodeGTF_Human protCodeGeneSym
elseif isequal(species{1},'Mouse')
    load GencodeGTF_Mouse protCodeGeneSym
end

pseudoInd_TypA=[];
for i=1:size(TypAExp_rmRep_dm,1)
    ind=find(strcmp(TypAExp_rmRep_dm.RowNames(i),protCodeGeneSym));
    if isempty(ind)
        pseudoInd_TypA=[pseudoInd_TypA;i];
    end
end
if ~isempty(pseudoInd_TypA)
    TypAExp_rmRep_dm(pseudoInd_TypA,:)=[]; %-checked.
end
TypAExp_rmRep_dmNew=TypAExp_rmRep_dm;
save Exp_cleaned_2 TypAExp_rmRep_dmNew -v7.3


%-3) Filter out lowly-expressed genes in cell type B.
%Data cleaning-1.
clear
TypBExp_rmRep_dm=bioma.data.DataMatrix('File', 'scRNAseq_CellTypeB.csv','Delimiter',',');
TypBExp_rmRep=double(TypBExp_rmRep_dm);
filterGene_TypB=[];
cutoff_TypB=floor(size(TypBExp_rmRep,2)*0.1);  %-must express in more than 10% (CellPhoneDB) of the cells in a specific cluster.
for p=1:size(TypBExp_rmRep,1)
    if length(find(TypBExp_rmRep(p,:)))<=cutoff_TypB
        filterGene_TypB=[filterGene_TypB;p];
    end
end
if ~isempty(filterGene_TypB)
    TypBExp_rmRep_dm(filterGene_TypB,:)=[]; %-checked.
end
save Exp_cleaned_3 TypBExp_rmRep_dm -v7.3

%-4) Only keep protein-coding genes.
%Data cleaning-2.
clear
load Exp_cleaned_3 TypBExp_rmRep_dm
fid0=fopen('Species.txt');
species=textscan(fid0,'%s','delimiter','\t');
fclose(fid0);
species=species{1};
if isequal(species{1},'Human')
    load GencodeGTF_Human protCodeGeneSym
elseif isequal(species{1},'Mouse')
    load GencodeGTF_Mouse protCodeGeneSym
end

pseudoInd_TypB=[];
for i=1:size(TypBExp_rmRep_dm,1)
    ind=find(strcmp(TypBExp_rmRep_dm.RowNames(i),protCodeGeneSym));
    if isempty(ind)
        pseudoInd_TypB=[pseudoInd_TypB;i];
    end
end
if ~isempty(pseudoInd_TypB)
    TypBExp_rmRep_dm(pseudoInd_TypB,:)=[]; %-checked.
end
TypBExp_rmRep_dmNew=TypBExp_rmRep_dm;
save Exp_cleaned_4 TypBExp_rmRep_dmNew -v7.3

%--Export gene expression into csv to compute
%mutual information-based indirect edge-filtered network matrices.
clear
load Exp_cleaned_2
load Exp_cleaned_4
emptyCharVector={''};
TypAExp_rmRep_dmNew=set(TypAExp_rmRep_dmNew,'Name',emptyCharVector); %formating the DataMatrix object.
TypBExp_rmRep_dmNew=set(TypBExp_rmRep_dmNew,'Name',emptyCharVector); 
dmwrite(TypAExp_rmRep_dmNew,'TypAExp_rmRepGf_dm.csv','Delimiter',',');
dmwrite(TypBExp_rmRep_dmNew,'TypBExp_rmRepGf_dm.csv','Delimiter',',');

exit  %-this is to exit MATLAB in order to run subsequent programs.


