function qualityList = gen_qualityCheck_v2(beta,omega,PCSF_nodeSym,PCSF_edgeSym,crosstalk_edgeSym,nTree,isolatedNode)
%GEN_QUALITYCHECK Summary of this function goes here
%   Detailed explanation goes here
qualityList=zeros(1,15); %-15 columns.

qualityList(1)=beta;
qualityList(2)=omega;

PCSF_geneSym=PCSF_nodeSym;
TypAGeneInd=[];
TypBGeneInd=[];
for i=1:length(PCSF_geneSym)
    
    if strcmp(PCSF_geneSym{i}(end-5:end),'__TypA') %-if true, it's TypA.
        TypAGeneInd=[TypAGeneInd;i];
        
    elseif strcmp(PCSF_geneSym{i}(end-5:end),'__TypB') %-if true, it's TypB.
        TypBGeneInd=[TypBGeneInd;i];
    end
end %-end of for

qualityList(3)=length(TypAGeneInd);
qualityList(4)=length(TypBGeneInd);
qualityList(5)=size(crosstalk_edgeSym,1);
qualityList(6)=size(crosstalk_edgeSym,1)/size(PCSF_edgeSym,1);

%---Fraction of crosstalk edges that have downstream intracellular
%pathways.
%1) construct adjacency matrix of PCSF.
adjMatrix=zeros(length(PCSF_geneSym),length(PCSF_geneSym));
for k=1:size(PCSF_edgeSym,1)
    ind_1=find(strcmp(PCSF_edgeSym(k,1),PCSF_geneSym));
    ind_2=find(strcmp(PCSF_edgeSym(k,2),PCSF_geneSym));
    adjMatrix(ind_1,ind_2)=1;
    adjMatrix(ind_2,ind_1)=1;
end
%2) check each crosstalk edge.
crossEdgeHasDown_TypA=0;
crossEdgeHasDown_TypB=0;
crossEdgeHasDown_TypATypB=0;
for s=1:size(crosstalk_edgeSym,1)
    ind_1=find(strcmp(crosstalk_edgeSym(s,1),PCSF_geneSym));
    ind_2=find(strcmp(crosstalk_edgeSym(s,2),PCSF_geneSym));
    flag_1=0;
    flag_2=0;
    
    %-ind_1 is TypA and ind_2 is TypB.
    if ismember(ind_1,TypAGeneInd) && ismember(ind_2,TypBGeneInd)
        ind_1_neighbor=find(adjMatrix(ind_1,:));
        if ~isempty(intersect(ind_1_neighbor,TypAGeneInd))
            crossEdgeHasDown_TypA=crossEdgeHasDown_TypA+1;
            flag_1=1;
        end
        
        ind_2_neighbor=find(adjMatrix(ind_2,:));
        if ~isempty(intersect(ind_2_neighbor,TypBGeneInd))
            crossEdgeHasDown_TypB=crossEdgeHasDown_TypB+1;
            flag_2=1;
        end
        
        if flag_1==1 && flag_2==1
            crossEdgeHasDown_TypATypB=crossEdgeHasDown_TypATypB+1;
        end
    end
   
    %-ind_1 is TypB and ind_2 is TypA.
    if ismember(ind_1,TypBGeneInd) && ismember(ind_2,TypAGeneInd)
        ind_1_neighbor=find(adjMatrix(ind_1,:));
        if ~isempty(intersect(ind_1_neighbor,TypBGeneInd))
            crossEdgeHasDown_TypB=crossEdgeHasDown_TypB+1;
            flag_1=1;
        end
        
        ind_2_neighbor=find(adjMatrix(ind_2,:));
        if ~isempty(intersect(ind_2_neighbor,TypAGeneInd))
            crossEdgeHasDown_TypA=crossEdgeHasDown_TypA+1;
            flag_2=1;
        end
        
        if flag_1==1 && flag_2==1
            crossEdgeHasDown_TypATypB=crossEdgeHasDown_TypATypB+1;
        end
    end
    
end %-end of for.
qualityList(7)=crossEdgeHasDown_TypA/size(crosstalk_edgeSym,1);
qualityList(8)=crossEdgeHasDown_TypB/size(crosstalk_edgeSym,1);
qualityList(9)=crossEdgeHasDown_TypATypB/size(crosstalk_edgeSym,1);

%---Number of hub nodes (degree > 10), trees and isolated nodes.
nodeDegree=sum(adjMatrix,2);
qualityList(10)=length(find(nodeDegree>10));
qualityList(11)=nTree; %-number of connected components.
qualityList(12)=length(isolatedNode);
qualityList(13)=nTree-length(isolatedNode); %-number of connected components with more than one node.

%-PCSF edge density to measure PCSF connectivity related to communication.
if length(PCSF_nodeSym)>=2
    possibleNum=nchoosek(length(PCSF_nodeSym),2);
    qualityList(14)=size(PCSF_edgeSym,1)/possibleNum;
else
    qualityList(14)=0;
end

%-Record PCSF edge number to determine min of omega.
qualityList(15)=size(PCSF_edgeSym,1); 

end


