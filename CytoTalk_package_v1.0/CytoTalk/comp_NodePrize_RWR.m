%------------Random walk with restart on TypA network.
%-Compute normalized adjacency matrix. 6min
clear
load MI_TypA MiMatrix_TypA
D_TypA=zeros(size(MiMatrix_TypA)); %-construct diagonal (weighted) degree matrix.
for k=1:size(MiMatrix_TypA,1)
    D_TypA(k,k)=sum(MiMatrix_TypA(k,:));
end
fix(clock)
Wnorm_TypA=(D_TypA^-0.5)*MiMatrix_TypA*(D_TypA^-0.5); %-checked. Result in symmetric matrix.
fix(clock)
save GeneRelevanceTypA Wnorm_TypA -v7.3

%-Extract all ligands or receptors as seeds for network propagation and construct prior information vector.
clear
load MI_TypA RowNames_TypA
fid0=fopen('Species.txt');
species=textscan(fid0,'%s','delimiter','\t');
fclose(fid0);
species=species{1};
if isequal(species{1},'Human')
    fid2=fopen('LigandReceptor_Human.txt');
elseif isequal(species{1},'Mouse')
    fid2=fopen('LigandReceptor_Mouse.txt');
end
LRsym=textscan(fid2,'%s %s','delimiter','\t');
fclose(fid2);
a=LRsym{1};
b=LRsym{2};
LR_sym=unique([a;b]);

Y_TypA=zeros(length(RowNames_TypA),1);
for i=1:length(RowNames_TypA)
    if ismember(RowNames_TypA(i),LR_sym)
        Y_TypA(i)=1;
    end
end
save GeneRelevanceTypA Y_TypA -append

%-Compute relevance value for each gene using random walk with restart procedure.
clear
load GeneRelevanceTypA Y_TypA Wnorm_TypA
iterNum=50; %enough to converge.
alpha=0.9;
Relevance_TypA=Y_TypA; %-initialization.
ReleDiff_TypA=ones(iterNum,1);
for i=1:iterNum
    %---------progress bar-------------%
    fprintf('Iteration %d.\n',i); %very fast.
    %----------------------------------%
    oldVec=Relevance_TypA;
    Relevance_TypA=alpha.*(Wnorm_TypA*Relevance_TypA)+(1-alpha).*Y_TypA; %-checked each step.
    
    %-compute difference between consecutive runs using 1-norm of a vector.
    ReleDiff_TypA(i)=norm(Relevance_TypA-oldVec,1)/norm(oldVec,1);
end
save GeneRelevanceTypA alpha Relevance_TypA ReleDiff_TypA -append


%------------Random walk with restart on TypB network.
%-Compute normalized adjacency matrix. 6min
clear
load MI_TypB MiMatrix_TypB
D_TypB=zeros(size(MiMatrix_TypB)); %-construct diagonal (weighted) degree matrix.
for k=1:size(MiMatrix_TypB,1)
    D_TypB(k,k)=sum(MiMatrix_TypB(k,:));
end
fix(clock)
Wnorm_TypB=(D_TypB^-0.5)*MiMatrix_TypB*(D_TypB^-0.5); %-checked. Result in symmetric matrix.
fix(clock)
save GeneRelevanceTypB Wnorm_TypB -v7.3

%-Extract all ligands or receptors as seeds for network propagation and construct prior information vector.
clear
load MI_TypB RowNames_TypB
fid0=fopen('Species.txt');
species=textscan(fid0,'%s','delimiter','\t');
fclose(fid0);
species=species{1};
if isequal(species{1},'Human')
    fid2=fopen('LigandReceptor_Human.txt');
elseif isequal(species{1},'Mouse')
    fid2=fopen('LigandReceptor_Mouse.txt');
end
LRsym=textscan(fid2,'%s %s','delimiter','\t');
fclose(fid2);
a=LRsym{1};
b=LRsym{2};
LR_sym=unique([a;b]);

Y_TypB=zeros(length(RowNames_TypB),1);
for i=1:length(RowNames_TypB)
    if ismember(RowNames_TypB(i),LR_sym)
        Y_TypB(i)=1;
    end
end
save GeneRelevanceTypB Y_TypB -append

%-Compute relevance value for each gene using random walk with restart procedure.
clear
load GeneRelevanceTypB Y_TypB Wnorm_TypB
iterNum=50;
alpha=0.9;
Relevance_TypB=Y_TypB; %-initialization.
ReleDiff_TypB=ones(iterNum,1);
for i=1:iterNum
    %---------progress bar-------------%
    fprintf('Iteration %d.\n',i); %very fast.
    %----------------------------------%
    oldVec=Relevance_TypB;
    Relevance_TypB=alpha.*(Wnorm_TypB*Relevance_TypB)+(1-alpha).*Y_TypB; %-checked each step.
    
    %-compute difference between consecutive runs using 1-norm of a vector.
    ReleDiff_TypB(i)=norm(Relevance_TypB-oldVec,1)/norm(oldVec,1);
end
save GeneRelevanceTypB alpha Relevance_TypB ReleDiff_TypB -append


