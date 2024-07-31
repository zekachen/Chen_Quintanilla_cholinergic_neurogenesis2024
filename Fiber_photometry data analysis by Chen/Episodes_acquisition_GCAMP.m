function [] = Episodes_acquisition_GCAMP

%% Episodes acquisition
% 
clc
all clear

Data=dir('10#REVERT_BASE_8HZ_p185_2_red__0__21-32-39-824.mat'); % extract the names of all txt files
path_data=pwd; % get current folder path
ii = 1:length(Data)
name= Data(ii).name (1:end-4);
load('10#REVERT_BASE_8HZ_p185_2_red__0__21-32-39-824.mat');
filter=0.02;

coG(1,:)=coef(1,5411:6610);
Gratiofilter(1,:)=filter_2sIIR(coG(1,:),filter,20,3,'high')+coG(1); 
coG(2,:)=coef(1,11731:12930);
Gratiofilter(2,:)=filter_2sIIR(coG(2,:),filter,20,3,'high')+coG(1); 
coG(3,:)=coef(1,18041:19240);
Gratiofilter(3,:)=filter_2sIIR(coG(3,:),filter,20,3,'high')+coG(1); 
coG(4,:)=coef(1,24361:25560);
Gratiofilter(4,:)=filter_2sIIR(coG(4,:),filter,20,3,'high')+coG(1); 
coG(5,:)=coef(1,30681:31880);
Gratiofilter(5,:)=filter_2sIIR(coG(5,:),filter,20,3,'high')+coG(1); 
coG(6,:)=coef(1,36991:38190);
Gratiofilter(6,:)=filter_2sIIR(coG(6,:),filter,20,3,'high')+coG(1); 

n=6
for i=1:n,
    meanCoGALLMEAN(i,1)=mean(Gratiofilter(i,1:1200));
    dFv2(i,:)=(Gratiofilter(i,:)-meanCoGALLMEAN(i,1))./meanCoGALLMEAN(i,1)*100;
    i=i+1;
 
end


%mean&AUC of dF/F_base_stim_post
n=6
for i=1:n,
    meandFbase(i,1)=mean(dFv2(i,1:300));
    meandFstim(i,1)=mean(dFv2(i,301:600));
    meandFpost(i,1)=mean(dFv2(i,601:900));

    AUCdFbase(i,1)=sum(dFv2(i,1:300))*0.1;
    AUCdFstim(i,1)=sum(dFv2(i,301:600))*0.1;
    AUCdFpost(i,1)=sum(dFv2(i,601:900))*0.1;
    
    i=i+1;
 
end

save(fullfile(path_data,[name '_AUCprocessed.mat']))
