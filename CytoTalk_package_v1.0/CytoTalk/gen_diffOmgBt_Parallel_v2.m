%--Generate beta and omega folders for PCSF running.
clear
load ArtiEdge artiEdge_id
load IntegratedNet_TypATypB_ID integratedNet_GenePrize_initial
fix(clock)
poolobj=parpool(4); %-this must be at most the real number of cores.
%-First beta range.
betaValues=5:5:300;  %-modify 1/4 for cancer study.
parfor idx=1:numel(betaValues)  %-modify 2/4 for cancer study.
    beta=betaValues(idx); %-modify 3/4 for cancer study.

    %-create "beta" outer folder.
    str1=strcat('./bt',num2str(beta,'%.6f'));
    status_1=mkdir(str1);
    
    %-generate prize txt file for this beta value.
    integratedNet_GenePrize=beta.*integratedNet_GenePrize_initial;
    
    PrizeFileName=strcat('./bt',num2str(beta,'%.6f'),'/IntegratedNet_nodePrize.txt');
    fid7=fopen(PrizeFileName,'w');
    formatSpec_2='%f\n';  %-the precision should be noted!!!
    for kk=1:length(integratedNet_GenePrize)
        fprintf(fid7,formatSpec_2,integratedNet_GenePrize(kk));
    end
    fclose(fid7);
    
    %-create "omega" inner folder.
    %-Second omega range.
    for omega=0.1:0.1:1.5 %-step by increment of middle.

        str2=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega));
        status_2=mkdir(str2);
        BtOmgFolderName=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega),'/');
        
        %--copy PCSF-input files.
        str3_source='IntegratedNet_edgeCost_common.txt';
        str3_dest=strcat(BtOmgFolderName,'IntegratedNet_edgeCost.txt');
        status_3=copyfile(str3_source,str3_dest);
        
        %--copy PCSF-generated files.
        str6_source='gen_PCSF.py';
        status_6=copyfile(str6_source,BtOmgFolderName);
        
        %--append different omega to the edge cost file in different folders.
        artiEdgeCost=omega*ones(size(artiEdge_id,1),1);
        EdgeCostFileName=strcat(BtOmgFolderName,'IntegratedNet_edgeCost.txt');
        fid20=fopen(EdgeCostFileName,'a+'); %-'a+' means append text to the file.
        formatSpec_20='%f\n';
        for kkk=1:size(artiEdgeCost,1)
            fprintf(fid20,formatSpec_20,artiEdgeCost(kkk));
        end
        fclose(fid20);
        
    end %-end of omega.
    
end %-end of beta.(parfor)
delete(poolobj);
fix(clock)


%--Generate job files for all PCSF running jobs.
clear
qqq=0; %-point to number of PCSF running jobs.
for beta=5:5:300  %-modify 4/4 for cancer study.

    for omega=0.1:0.1:1.5 %-step by increment of middle.
        
        BtOmgFolderName=strcat('./bt',num2str(beta,'%.6f'),'/bt',num2str(beta,'%.6f'),'_omg',num2str(omega),'/');

        qqq=qqq+1;
        cmd_1=strcat('cd',{' '},BtOmgFolderName); %-strcat will remove tailing blank space, but this case can reserve blank space.
        cmd_2='python gen_PCSF.py';
        JobFileName=strcat('PCSFjob_',num2str(qqq),'.sh');
        fid25=fopen(JobFileName,'w');
        
        fprintf(fid25,'#!/bin/bash'); %-first line.
        fprintf(fid25,'\n');
        fprintf(fid25,'\n'); %-second empty line.
        fprintf(fid25,cmd_1{1}); %-third line.
        fprintf(fid25,'\n');
        fprintf(fid25,cmd_2); %-fourth line.
        
        fclose(fid25);
        
    end %-end of omega.
    
end %-end of beta.

%--Generate a job summary file.
clear
JobTable=dir('PCSFjob_*.sh');
nJob=size(JobTable,1);
fid26=fopen('Jobs2Run','w');
for rrr=1:nJob
    cmd=strcat('bash PCSFjob_',num2str(rrr),'.sh');
    if rrr~=nJob
        fprintf(fid26,cmd);
        fprintf(fid26,'\n');
    else
        fprintf(fid26,cmd); %-final line doesn't have new line character.
    end
end
fclose(fid26);

%--Generate a bash file to parallel run PCSF jobs.
clear
fid28=fopen('JobRun_Parallel.sh','w');
        
fprintf(fid28,'#!/bin/bash'); %-first line.
fprintf(fid28,'\n');
fprintf(fid28,'\n'); %-second empty line.
fprintf(fid28,'date +%%D-%%H:%%M:%%S'); %-third line to show current time.
fprintf(fid28,'\n');
fprintf(fid28,'parallel --jobs 6 < Jobs2Run'); %-fourth line.
fprintf(fid28,'\n');
fprintf(fid28,'date +%%D-%%H:%%M:%%S'); %-fifth line to show current time.

fclose(fid28);


