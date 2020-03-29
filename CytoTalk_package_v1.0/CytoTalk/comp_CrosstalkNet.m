clear
load Crosstalk_TypATypB crosstalk_Sym crosstalk_Sco
load MI_topNet_TypATypB MiList_geneSym_TypATypB

crosstalkNet_Sym=cell(size(crosstalk_Sym,1),2); %-initiated as maximum.
crosstalkNet_id=zeros(size(crosstalk_Sym,1),2);
crosstalkNet_Sco=ones(size(crosstalk_Sym,1),1)*(-100); %-initiated as -100.

k=1; %-point to crosstalkNet_Sym.
for i=1:size(crosstalk_Sym,1)
    %---------progress bar-------------%
    %fprintf('crosstalk_Sym %d.\n',i);
    %----------------------------------%
    if ismember(crosstalk_Sym{i,1},MiList_geneSym_TypATypB) && ismember(crosstalk_Sym{i,2},MiList_geneSym_TypATypB)
        crosstalkNet_Sym{k,1}=crosstalk_Sym{i,1};
        crosstalkNet_Sym{k,2}=crosstalk_Sym{i,2};
            
        crosstalkNet_id(k,1)=find(strcmp(crosstalkNet_Sym{k,1},MiList_geneSym_TypATypB));
        crosstalkNet_id(k,2)=find(strcmp(crosstalkNet_Sym{k,2},MiList_geneSym_TypATypB));
            
        crosstalkNet_Sco(k)=crosstalk_Sco(i);

        k=k+1;
    end
end

%-clean "crosstalkNet_Sym" and "crosstalkNet_Sco".
ind=find(cellfun(@isempty,crosstalkNet_Sym(:,1)));
indind=find(crosstalkNet_Sco==-100);
if isequal(ind,indind)  %if true, correct!
    crosstalkNet_Sym(indind,:)=[];
    crosstalkNet_Sco(indind,:)=[];
    crosstalkNet_id(indind,:)=[];
end

%--Remove zero crosstalk edges.
negInd=find(crosstalkNet_Sco<=0);
crosstalkNet_Sym(negInd,:)=[];
crosstalkNet_Sco(negInd,:)=[];
crosstalkNet_id(negInd,:)=[];

%-edge weight normalized using zscore.
mu=mean(crosstalkNet_Sco);
sigma=std(crosstalkNet_Sco);
crosstalkNet_Sco=(crosstalkNet_Sco-mu)/sigma;
save CrosstalkNet_TypATypB crosstalkNet_Sym crosstalkNet_Sco crosstalkNet_id


