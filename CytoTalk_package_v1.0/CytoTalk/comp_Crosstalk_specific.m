%-------------------------Compute crosstalk score-----------------------%
%-import non-self-talk score in TypB.
clear
fid2=fopen('NonSelfTalkSym_TypB.txt');
LRsym=textscan(fid2,'%s %s','delimiter','\t');
fclose(fid2);
a=LRsym{1};
b=LRsym{2};
LRpairSym_TypB=[a,b];

fid3=fopen('NonSelfTalkSco_TypB.txt');
NonSelfTalkSco_TypB=textscan(fid3,'%f','TreatAsEmpty','NA');
fclose(fid3);
NonSelfTalkSco_TypB=NonSelfTalkSco_TypB{1};  
negativeInd=find(NonSelfTalkSco_TypB<0);
NonSelfTalkSco_TypB(negativeInd)=0; 
nanInd=find(isnan(NonSelfTalkSco_TypB));
NonSelfTalkSco_TypB(nanInd)=0;
save Crosstalk_TypATypB LRpairSym_TypB NonSelfTalkSco_TypB

%-import non-self-talk score in TypA.
clear
fid2=fopen('NonSelfTalkSym_TypA.txt');
LRsym=textscan(fid2,'%s %s','delimiter','\t');
fclose(fid2);
a=LRsym{1};
b=LRsym{2};
LRpairSym_TypA=[a,b]; 

fid3=fopen('NonSelfTalkSco_TypA.txt');
NonSelfTalkSco_TypA=textscan(fid3,'%f','TreatAsEmpty','NA');
fclose(fid3);
NonSelfTalkSco_TypA=NonSelfTalkSco_TypA{1};  
negativeInd=find(NonSelfTalkSco_TypA<0);
NonSelfTalkSco_TypA(negativeInd)=0; 
nanInd=find(isnan(NonSelfTalkSco_TypA));
NonSelfTalkSco_TypA(nanInd)=0;
save Crosstalk_TypATypB LRpairSym_TypA NonSelfTalkSco_TypA -append

clear
load Crosstalk_TypATypB
load GeneCellTypeSpecific
%TypAExp_rmRep_dm=bioma.data.DataMatrix('File', 'scRNAseq_CellTypeA.csv','Delimiter',',');
%TypBExp_rmRep_dm=bioma.data.DataMatrix('File', 'scRNAseq_CellTypeB.csv','Delimiter',',');

%TypB_avg=mean(double(TypBExp_rmRep_dm),2);
%TypA_avg=mean(double(TypAExp_rmRep_dm),2);

if isequal(LRpairSym_TypB,LRpairSym_TypA)  %-should same, otherwise, need preprocessing.
    LRpairSym=LRpairSym_TypA;
    crosstalk_Sym=cell(size(LRpairSym,1)*2,2); %-initiated as double edges (maximum).
    ExpressedSco=ones(size(LRpairSym,1)*2,1)*(-100); %-initiated as -100.
    NonSelfTalkSco=ones(size(LRpairSym,1)*2,1)*(-100); %-initiated as -100.
    k=1; %-point to crosstalk_Sym.
    for i=1:size(LRpairSym,1)
        %---------progress bar-------------%
        %fprintf('LRpairSym %d.\n',i);
        %----------------------------------%
        tf=strcmp(LRpairSym{i,1},LRpairSym{i,2});
        if ~tf
                %--From TypA to TypB.
                crosstalk_Sym{k,1}=strcat(LRpairSym{i,1},'__TypA'); %gene in TypA cell.
                crosstalk_Sym{k,2}=strcat(LRpairSym{i,2},'__TypB'); %gene in TypB cell.
                
                %-compute non-self-talk score.
                NonSelfTalkSco(k)=(NonSelfTalkSco_TypA(i)+NonSelfTalkSco_TypB(i))/2;
                
                %-compute expressed score.
                %ind_TypA=find(strcmp(LRpairSym{i,1},TypAExp_rmRep_dm.RowNames)); %-must have a single value.
                %exp_TypA=TypA_avg(ind_TypA);
                ind_TypA=find(strcmp(LRpairSym{i,1},same_rowname)); %-must have a single value.
                exp_TypA=typeSpecific{1}(ind_TypA);
                if exp_TypA>0
                    exp_TypA=exp_TypA; %no change.
                else
                    exp_TypA=0;
                end
                
                %ind_TypB=find(strcmp(LRpairSym{i,2},TypBExp_rmRep_dm.RowNames)); %-must have a single value.
                %exp_TypB=TypB_avg(ind_TypB);
                ind_TypB=find(strcmp(LRpairSym{i,2},same_rowname)); %-must have a single value.
                exp_TypB=typeSpecific{2}(ind_TypB);
                if exp_TypB>0
                    exp_TypB=exp_TypB; %no change.
                else
                    exp_TypB=0;
                end
 
                ExpressedSco(k)=(exp_TypA+exp_TypB)/2;

                k=k+1;
                
                %--From TypB to TypA.
                crosstalk_Sym{k,1}=strcat(LRpairSym{i,1},'__TypB'); %gene in TypB cell.
                crosstalk_Sym{k,2}=strcat(LRpairSym{i,2},'__TypA'); %gene in TypA cell.
                
                %-compute non-self-talk score.
                NonSelfTalkSco(k)=(NonSelfTalkSco_TypB(i)+NonSelfTalkSco_TypA(i))/2; %-same as before.
                
                %-compute expressed score.
                %ind_TypB=find(strcmp(LRpairSym{i,1},TypBExp_rmRep_dm.RowNames)); %-must have a single value.
                %exp_TypB=TypB_avg(ind_TypB);
                ind_TypB=find(strcmp(LRpairSym{i,1},same_rowname)); %-must have a single value.
                exp_TypB=typeSpecific{2}(ind_TypB);
                if exp_TypB>0
                    exp_TypB=exp_TypB; %no change.
                else
                    exp_TypB=0;
                end
                
                %ind_TypA=find(strcmp(LRpairSym{i,2},TypAExp_rmRep_dm.RowNames)); %-must have a single value.
                %exp_TypA=TypA_avg(ind_TypA);
                ind_TypA=find(strcmp(LRpairSym{i,2},same_rowname)); %-must have a single value.
                exp_TypA=typeSpecific{1}(ind_TypA);
                if exp_TypA>0
                    exp_TypA=exp_TypA; %no change.
                else
                    exp_TypA=0;
                end
                
                ExpressedSco(k)=(exp_TypB+exp_TypA)/2;

                k=k+1;
        
       else    %-two gene symbols are the same. No need to switch.
                %--From TypA to TypB.
                crosstalk_Sym{k,1}=strcat(LRpairSym{i,1},'__TypA'); %gene in TypA cell.
                crosstalk_Sym{k,2}=strcat(LRpairSym{i,2},'__TypB'); %gene in TypB cell.
                
                %-compute non-self-talk score.
                NonSelfTalkSco(k)=(NonSelfTalkSco_TypA(i)+NonSelfTalkSco_TypB(i))/2;
                
                %-compute expressed score.
                %ind_TypA=find(strcmp(LRpairSym{i,1},TypAExp_rmRep_dm.RowNames)); %-must have a single value.
                %exp_TypA=TypA_avg(ind_TypA);
                ind_TypA=find(strcmp(LRpairSym{i,1},same_rowname)); %-must have a single value.
                exp_TypA=typeSpecific{1}(ind_TypA);
                if exp_TypA>0
                    exp_TypA=exp_TypA; %no change.
                else
                    exp_TypA=0;
                end
                
                %ind_TypB=find(strcmp(LRpairSym{i,2},TypBExp_rmRep_dm.RowNames)); %-must have a single value.
                %exp_TypB=TypB_avg(ind_TypB);
                ind_TypB=find(strcmp(LRpairSym{i,2},same_rowname)); %-must have a single value.
                exp_TypB=typeSpecific{2}(ind_TypB);
                if exp_TypB>0
                    exp_TypB=exp_TypB; %no change.
                else
                    exp_TypB=0;
                end
                
                ExpressedSco(k)=(exp_TypA+exp_TypB)/2;

                k=k+1;

        end  %-end of if ~tf.
    end %-end of for.
    
end  %-end of if.

%-clean "crosstalk_Sym", "ExpressedSco" and "NonSelfTalkSco".
ind=find(cellfun(@isempty,crosstalk_Sym(:,1)));
indind=find(NonSelfTalkSco==-100);
if isequal(ind,indind)  %if true, correct!
    crosstalk_Sym(indind,:)=[];
    ExpressedSco(indind,:)=[];
    NonSelfTalkSco(indind,:)=[];
end
save Crosstalk_TypATypB crosstalk_Sym ExpressedSco NonSelfTalkSco -append

%--Compute crosstalk score using normalized ExpressedSco and Non-self-talk score.
clear
load Crosstalk_TypATypB ExpressedSco NonSelfTalkSco
%-normalize ExpressedSco.
ExpressedSco_norm=(ExpressedSco-min(ExpressedSco))/(max(ExpressedSco)-min(ExpressedSco));

%-normalize NonSelfTalkSco.
NonSelfTalkSco_norm=(NonSelfTalkSco-min(NonSelfTalkSco))/(max(NonSelfTalkSco)-min(NonSelfTalkSco));

%-compute crosstalk score.
crosstalk_Sco=ExpressedSco_norm.*NonSelfTalkSco_norm;
save Crosstalk_TypATypB crosstalk_Sco -append


