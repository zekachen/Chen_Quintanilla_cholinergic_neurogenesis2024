function [] = unmixing_nofilterAUC_zeka2023
%%
% Read raw mixed spectra time series from Oceanview output, and unmix the spectra time series one by one using user defined spectra reference.
%%

clc
all clear
Data=dir('*.txt'); % extract the names of all txt files
path_data=pwd; % get current folder path
[refID,path_ref] = uigetfile('*.csv','Select reference');
ref=csvread([path_ref refID],1,1);
for ii = 1:length(Data)
    i=0; test = {{'a'}}; % skipping headers
    while isnan(str2double(test{1,1})) == 1
        file = fopen([path_data,'/',Data(ii).name],'r');
        name = Data(ii).name (1:end-4);  
        test = textscan(file, '%s',1,'HeaderLines',i);
        fclose(file);
        i=i+1;
    end
    if i == 2
        i = i+1; 
    else
        i = i;
    end
    class=[];
    for j=1:1044
        class=[class '%f '];
    end

    file = fopen([path_data,'/',Data(ii).name],'r');
    data = textscan(file, ['%s' '%s' '%s' class],'HeaderLines',i); 
    data = cell2mat(data(4:end));  % remove the first 3 columns
    data = data';
    fclose(file);
    coef=zeros(size(ref,2),size(data,2));
    for j=1:size(data,2)
        coef(:,j)=max(0,lsqnonneg(ref(140:500,:), data(140:500,j)));
        clc
        [num2str(j/size(data,2)*100) '%']
    end
    figure
    for j=1:size(coef,1)
        subplot(size(coef,1)+1,1,j)
        plot(0.1:0.1:size(data,2)/10,coef(j,:)/mean(coef(j,:),2)*100-100)
        xlabel('Time (s)','FontWeight','bold','FontSize',12)
        ylabel('dF/F (%)','FontWeight','bold','FontSize',12)
    end
    subplot(size(coef,1)+1,1,j+1)
    plot(0.1:0.1:size(data,2)/10,coef(1,:)./coef(2,:))
    title('y=coG coTd ratio')
        xlabel('Time (s)','FontWeight','bold','FontSize',12)
        ylabel('Ratio','FontWeight','bold','FontSize',12)
    save([Data(ii).name(1:length(Data(ii).name)-4) '.mat'],'coef')
    savefig([Data(ii).name(1:length(Data(ii).name)-4) '.fig'])
    exportgraphics(gcf, [Data(ii).name(1:length(Data(ii).name)-4) '.jpg']);
end

%%  definate fs(frequency), sd(Standard deviation)  
prompt = {'fs','sd'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'10','2'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
fs=str2num(answer{1});  
sd=str2num(answer{2});

%%  no_filter_for AUC analysis
prompt = {'Enter begining_base(10*s):','Enter finishing_base(10*s):','Enter begining_event(10*s):','Enter finishing_event(10*s):','Enter begining_after(10*s):','Enter finishing_after(10*s):'};%select episodes, make sure the length of each episode is same.
dlgtitle = 'Input';
dims = [1 35];
definput = {'500','2500','3500','5500','6500','8500'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
T1=str2num(answer{1});
T2=str2num(answer{2});
T3=str2num(answer{3});
T4=str2num(answer{4});
T5=str2num(answer{5});
T6=str2num(answer{6});
coG=coef(1,1:8979); %the time length based on your recording 
coG_base=coG(1,T1:1:T2);
mean_coG_base=mean(coG_base);
coGratio=(coG-mean_coG_base)./mean_coG_base*100;%to get deltaF/F

%%
raw_base=coGratio(1,T1:1:T2);
mean_raw_base=mean(raw_base);
raw_event=coGratio(1,T3:1:T4);
mean_raw_event=mean(raw_event);
raw_after=coGratio(1,T5:1:T6);
mean_raw_after=mean(raw_after);
mean_sum=mean(coGratio)
s = size(coGratio,2);   % length of the recording
time = [1:1:s]./60./fs;  % convert to minute

%% base threshold analysis
% thresholding total recording session
% find peak raw data, thresholding at 2SD
     [pks_base,locs_base] = findpeaks(raw_base,'MinPeakDistance',1);
      threshold_base_SD = sd*std(raw_base);
      locs_base_SD = locs_base(pks_base>threshold_base_SD);
      pks_base_SD = raw_base(locs_base_SD);

% exclude peaks from 1st round, and find peak for rest of the data, thresholding at 2SD
    base_ratio_2 = raw_base;
    base_ratio_2 (base_ratio_2 > threshold_base_SD) = nan;
    threshold_base_SD_2 = sd*std(base_ratio_2(~isnan(base_ratio_2)));
    [pks_base_2,locs_base_2] = findpeaks(base_ratio_2,'MinPeakDistance',1);
    locs_base_SD_2 = locs_base_2(pks_base_2>threshold_base_SD_2);
    pks_base_SD_2 = raw_base(locs_base_SD_2);

% exclude peaks from 1st, 2nd round, and find peak for rest of the data, thresholding at 2SD
    base_ratio_3 = raw_base;
    base_ratio_3 (base_ratio_3 > threshold_base_SD_2) = nan;
    threshold_base_SD_3 = sd*std(base_ratio_3(~isnan(base_ratio_3)));
    [pks_base_3,locs_base_3] = findpeaks(base_ratio_3,'MinPeakDistance',1);
    locs_base_SD_3 = locs_base_3(pks_base_3>threshold_base_SD_3);
    pks_base_SD_3 = raw_base(locs_base_SD_3);

AUC_base = sum(raw_base(raw_base>threshold_base_SD_3)); % sum of everything above the 3rd threshold
above_base_base = sum(raw_base>threshold_base_SD_3)./s*100; % sum of datapoints above the 3rd threshold

AUC_event = sum(raw_event(raw_event>threshold_base_SD_3)); % sum of everything above the 3rd threshold
above_base_event = sum(raw_event>threshold_base_SD_3)./s*100; % sum of datapoints above the 3rd threshold

AUC_after = sum(raw_after(raw_after>threshold_base_SD_3)); % sum of everything above the 3rd threshold
above_base_after = sum(raw_after>threshold_base_SD_3)./s*100; % sum of datapoints above the 3rd threshold


%% plot calcium traces
    figure
    x0=100;
    y0=200;
    width=1200;
    height=500;
    set(gcf,'position',[x0,y0,width,height])   % set figure configurations    

%% Base_2SD_above
    subplot(2,1,1)
% plot traces and detected calcium events
plot (time,coGratio, 'k')
    hold on
% plot(coGratio,'k')
threshold_base_SD_seq_3 = ones(1, s).*threshold_base_SD_3;
plot(time,threshold_base_SD_seq_3,'LineWidth',2, 'Color', [0, 0, 1, 0.2])

 title('AUC above 2SD')
 xlabel('Time (mins)')
 ylabel('dF/F%')

%% Base_mean(base)_above
    subplot(2,1,2)
    plot (time,coGratio, 'k')
    hold on
%     raw_base_line = sd*std(raw_base);
    threshold_mean_base = ones(1, s).*mean_raw_base;
    plot(time,threshold_mean_base,'LineWidth',2, 'Color', [1, 0, 0, 0.2])

    title('AUC above mean')
    xlabel('Time (mins)')
    ylabel('dF/F%')
%     
%% save data and variables to the summarized file
    savefig([Data(ii).name(1:length(Data(ii).name)-4) '_AUC_event.fig'])
    saveas(gcf,fullfile(path_data,[name '_AUC_event.png']))
    exportgraphics(gcf,fullfile(path_data,[name '_AUC_event.pdf']),'ContentType','vector')

save(fullfile(path_data,[name '_AUCprocessed.mat']))

   
end