%-Generate integrated network for PCSF.
clear
load MI_topNet_TypATypB
load CrosstalkNet_TypATypB
load ArtiEdge
integratedNet_EdgeID=[MiList_genePair_TypATypB;crosstalkNet_id;artiEdge_id];
integratedNet_EdgeID=integratedNet_EdgeID-1; %-Python is 0-based indexing.
integratedNet_EdgeValue_common=[MiList_value_TypATypB;crosstalkNet_Sco]; %-need add artiEdgeCost.

minValue=min(integratedNet_EdgeValue_common);
maxValue=max(integratedNet_EdgeValue_common);
integratedNet_EdgeValue_common=(integratedNet_EdgeValue_common-minValue)/(maxValue-minValue);
integratedNet_EdgeCost_common=1-integratedNet_EdgeValue_common;

integratedNet_GeneSym=[MiList_geneSym_TypATypB;artiNodeSym];
integratedNet_GenePrize_initial=[MiList_genePrize_TypATypB;artiNodePrize]; %-need multiply beta.
save IntegratedNet_TypATypB_ID integratedNet_EdgeID integratedNet_EdgeCost_common integratedNet_GeneSym integratedNet_GenePrize_initial


