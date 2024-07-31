function [] = unmixing_WV456_540_2023zekaedit
%%
% Read raw mixed spectra time series from Oceanview output, and unmix the spectra time series one by one using user defined spectra reference.
% Tzu-Hao Harry Chao 2020/02/14
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
%         file = fopen([path_data dataID],'r');
        test = textscan(file, '%s',1,'HeaderLines',i);
        fclose(file);
        i=i+1;
    end
    
    % to fix additional gap line before the wavelength column
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
        coef(:,j)=max(0,lsqnonneg(ref(140:250,:), data(140:250,j)));% acquisition of siganls of 456_540 nm wavelenth 
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
