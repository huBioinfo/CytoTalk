%-Add all gene-related artificial edges.
clear
load MI_topNet_TypATypB
artiConnectGene_id=unique([MiList_genePair_TypATypB(:,1);MiList_genePair_TypATypB(:,2)]);

artiNodePrize=1; %-This value is irrelavant to the result.
artiNodeSym={'artiNode'};
artiNode_id=length(MiList_geneSym_TypATypB)+1;
artiConnectGene_id=unique(artiConnectGene_id);
artiEdge_id=[artiNode_id.*ones(length(artiConnectGene_id),1),artiConnectGene_id];
save ArtiEdge artiNodePrize artiNodeSym artiNode_id artiEdge_id


