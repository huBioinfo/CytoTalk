%---Import edge indices and convert them to edge symbols.
clearvars -except beta omega
fid4=fopen('PCSF_Edge.txt');
PCSF_edgeIndex=textscan(fid4,'%d');
PCSF_edgeIndex=PCSF_edgeIndex{1}; %-checked.
fclose(fid4);

fid5=fopen('PCSF_Node.txt');
PCSF_nodeIndex=textscan(fid5,'%d');
PCSF_nodeIndex=PCSF_nodeIndex{1}; %-checked.
fclose(fid5);
save pcsf_result4 PCSF_edgeIndex PCSF_nodeIndex

clearvars -except beta omega
load pcsf_result4
load ../../IntegratedNet_TypATypB_ID
PCSF_edgeSym=cell(length(PCSF_edgeIndex),2);
for i=1:length(PCSF_edgeIndex)
    PCSF_edgeSym{i,1}=integratedNet_GeneSym{integratedNet_EdgeID(PCSF_edgeIndex(i)+1,1)+1}; %-Note:edgeIndex and nodeIndex both from 0 in Python.
    PCSF_edgeSym{i,2}=integratedNet_GeneSym{integratedNet_EdgeID(PCSF_edgeIndex(i)+1,2)+1};
end
PCSF_edge2nodeSym=unique([PCSF_edgeSym(:,1);PCSF_edgeSym(:,2)]); %-this is still a tree.
PCSF_nodeSym=integratedNet_GeneSym(PCSF_nodeIndex+1);
if isequal(unique(PCSF_edge2nodeSym),unique(PCSF_nodeSym)) %-"unique" can sort.
    correctFlag=1;
else
    correctFlag=0;
end
save pcsf_result4 PCSF_edgeSym PCSF_nodeSym correctFlag -append

%---Remove 'artiNode'-associated edges.
clearvars -except beta omega
load pcsf_result4
load ../../IntegratedNet_TypATypB_ID
removeInd=[];
Root=[];
for i=1:size(PCSF_edgeSym,1)
    ind=find(strcmp('artiNode',PCSF_edgeSym(i,:)));
    if ~isempty(ind)
        removeInd=[removeInd;i];
        
        if ind==1  %-record root nodes.
            Root=[Root;PCSF_edgeSym(i,2)];
        elseif ind==2
            Root=[Root;PCSF_edgeSym(i,1)];
        end
        
    end
end
PCSF_edgeIndex(removeInd,:)=[]; %-exclude artificial edge index.
PCSF_edgeSym(removeInd,:)=[];
PCSF_edgeCost=integratedNet_EdgeCost_common(PCSF_edgeIndex+1); %-this already exclude artifical edge costs.
nEdge=size(PCSF_edgeSym,1);
nTree=length(Root);
save pcsf_result4 PCSF_edgeIndex PCSF_edgeSym PCSF_edgeCost Root nEdge nTree -append

%---Remove 'artiNode' from PCSF_nodeSym and then find isolated nodes.
clearvars -except beta omega
load pcsf_result4
ind=find(strcmp('artiNode',PCSF_nodeSym));
PCSF_nodeSym(ind,:)=[];
PCSF_nodeIndex(ind,:)=[];
isolatedNode=setdiff(PCSF_nodeSym,unique(PCSF_edgeSym(:)));
save pcsf_result4 PCSF_nodeSym PCSF_nodeIndex isolatedNode -append


%---Extract # of genes in different cell types.
clearvars -except beta omega
load pcsf_result4 PCSF_nodeSym
PCSF_geneSym=PCSF_nodeSym;
nTypAGene=0;
nTypBGene=0;
for i=1:length(PCSF_geneSym)
    %---------progress bar-------------%
%     fprintf('PCSF_geneSym %d.\n',i);
    %----------------------------------%
    if strcmp(PCSF_geneSym{i}(end-5:end),'__TypA') %-if true, it's TypA.
        nTypAGene=nTypAGene+1;
        
    elseif strcmp(PCSF_geneSym{i}(end-5:end),'__TypB') %-if true, it's TypB.
        nTypBGene=nTypBGene+1;
    end
    
end
save pcsf_result4 nTypAGene nTypBGene -append


%-----------Extract crosstalk edges in the forest.
clearvars -except beta omega
load pcsf_result4 PCSF_edgeSym
PCSF_edgeCellType=zeros(size(PCSF_edgeSym,1),1);
for i=1:size(PCSF_edgeSym,1)
    if strcmp(PCSF_edgeSym{i,1}(end-5:end),'__TypA') && strcmp(PCSF_edgeSym{i,2}(end-5:end),'__TypA')
        PCSF_edgeCellType(i)=1; %-"1" as TypA.
    elseif strcmp(PCSF_edgeSym{i,1}(end-5:end),'__TypB') && strcmp(PCSF_edgeSym{i,2}(end-5:end),'__TypB')
        PCSF_edgeCellType(i)=2; %-"2" as TypB.
    else %-must be crosstalk type.
        PCSF_edgeCellType(i)=3; %-"3" as Crosstalk.
    end
end
crosstalkInd=find(PCSF_edgeCellType==3);
crosstalk_edgeSym=PCSF_edgeSym(crosstalkInd,:);
crosstalk_edgeNum=size(crosstalk_edgeSym,1);
save pcsf_result4 crosstalk_edgeSym crosstalk_edgeNum PCSF_edgeCellType -append


%-Extract non-crosstalk LR pairs in the network.
clearvars -except beta omega
load pcsf_result4 crosstalk_edgeSym
load ../../CrosstalkNet_TypATypB crosstalkNet_Sym crosstalkNet_Sco
crosstalk_edgeInd=zeros(size(crosstalk_edgeSym,1),1);
wrongFlag2_crosstalk=0;
wrongCrosstalk2_edgeInd=[];
for i=1:size(crosstalk_edgeSym,1)
    %---------progress bar-------------%
%     fprintf('crosstalk_edgeSym %d.\n',i);
    %----------------------------------%
    vec_11=strcmp(crosstalk_edgeSym{i,1},crosstalkNet_Sym(:,1));
    vec_22=strcmp(crosstalk_edgeSym{i,2},crosstalkNet_Sym(:,2));
    ind_a=find(vec_11 & vec_22); %-must single value or empty.
    
    vec_12=strcmp(crosstalk_edgeSym{i,1},crosstalkNet_Sym(:,2));
    vec_21=strcmp(crosstalk_edgeSym{i,2},crosstalkNet_Sym(:,1));
    ind_b=find(vec_12 & vec_21); %-must single value or empty.
    
    if ~xor(isempty(ind_a),isempty(ind_b)) %if both true or both false, then set flag.
        wrongFlag2_crosstalk=1;
        wrongCrosstalk2_edgeInd=[wrongCrosstalk2_edgeInd;i];
    end
    
    if ~isempty(ind_a)
        crosstalk_edgeInd(i)=ind_a; %-checked.
    elseif ~isempty(ind_b)
        crosstalk_edgeInd(i)=ind_b;
    end
end
noncrosstalk_edgeInd=setdiff([1:size(crosstalkNet_Sym,1)]',crosstalk_edgeInd);
noncrosstalk_edgeSym=crosstalkNet_Sym(noncrosstalk_edgeInd,:);
noncrosstalk_edgeNum=size(noncrosstalk_edgeSym,1);

%-Extracting ORIGINAL crosstalk score of crosstalk edges and non-crosstalk edges.
crosstalk_edgeSco=crosstalkNet_Sco(crosstalk_edgeInd);
noncrosstalk_edgeSco=crosstalkNet_Sco(noncrosstalk_edgeInd);
save pcsf_result4 noncrosstalk_edgeSym wrongFlag2_crosstalk wrongCrosstalk2_edgeInd noncrosstalk_edgeNum crosstalk_edgeSco noncrosstalk_edgeSco -append

%-Compute the one-sided rank-sum test between crosstak and non-crosstalk LR
%pairs.
clearvars -except beta omega
load pcsf_result4 crosstalk_edgeSco noncrosstalk_edgeSco
if ~isempty(crosstalk_edgeSco) && ~isempty(noncrosstalk_edgeSco)
    p_crosstalkSco=ranksum(crosstalk_edgeSco,noncrosstalk_edgeSco,'tail','right');
else
    p_crosstalkSco=1;
end
save pcsf_result4 p_crosstalkSco -append


