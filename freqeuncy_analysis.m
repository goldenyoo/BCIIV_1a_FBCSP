%% Call raw data
close all
clear all

% Ask user for input parameters
prompt = {'Data label: ', 'Feature vector length: ', 'Re-referencing: 0 (Non),1 (CAR), 2 (LAP)', 'BFP order'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'a', '3', '1','20'};
answer = inputdlg(prompt,dlgtitle,dims,definput);


% Error detection
if isempty(answer), error("Not enough input parameters."); end

% Input parameters
data_label = string(answer(1,1));   % Calib_ds1 + "data_label"
m = double(string(answer(2,1))); % feature vector will have length (2m)
referencing = double(string(answer(3,1)));
order = double(string(answer(4,1)));

tmp = ["a", "b", "f", "g"];

% for label = 1: length(tmp)
%     data_label = tmp(label);

% Load file
FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_calib_ds1',data_label,'.mat');
load(FILENAME);

% Data rescale
cnt= 0.1*double(cnt);
cnt = cnt';

% Use designated electrode (C3, C4)
cnt_c = cnt(3:55,:);

%% Preprocessing
if referencing ~= 0
    %%% Calculate differential voltage
    for i = 1 : size(cnt_c,1)
        cnt_c(i,:) = cnt_c(i,:) - cnt(29,:);
    end

    
    if referencing == 1 % common average
        Means = (1/size(cnt_c,1))*sum(cnt_c);
        for i = 1 : size(cnt_c,1)
            cnt_c(i,:) = cnt_c(i,:) - Means;
        end
    elseif referencing == 2 % LAP
        %%%
    end
end



%% 
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
% for i = 1:length(mrk.pos)  
for i = 1:4 
    
    % One trial data
    E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+350); 
    E = E.^2;
    E = sum(E);
    
    %Plot freqeuncy analysis
    if i == 1
        figure(f1);        
     elseif i==2
        figure(f2);
     elseif i==3
         figure(f3);
     else
         figure(f4);
    end 
    subplot(2,1,1);
    sinad(E(1,:),100);
%     subplot(2,2,3);
%     sinad(E(2,:),100);
    % According to its class, divide calculated covariance
    if mrk.y(1,i) == 1 
        %
    else
        %
    end
end


%% 
%BPF Design
bpFilt = designfilt('bandpassfir','FilterOrder',order, ...
         'CutoffFrequency1',8,'CutoffFrequency2',30, ...
         'SampleRate',100);

% Apply BPF
for i = 1:size(cnt_c,1)
    cnt_c(i,:) = filtfilt(bpFilt, cnt_c(i,:));
end
%% 
% for i = 1:length(mrk.pos)  
for i = 1:4 
    
    % One trial data
    E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+350);   
    E = E.^2;
    E = sum(E);
    
    %Plot freqeuncy analysis
     if i == 1
        figure(f1);        
     elseif i==2
        figure(f2);
     elseif i==3
         figure(f3);
     else
         figure(f4);
    end 
    subplot(2,1,2);
    sinad(E(1,:),100);
%     subplot(2,2,4);
%     sinad(E(2,:),100);
    % According to its class, divide calculated covariance
    if mrk.y(1,i) == 1 
        %
    else
        %
    end
end
% end
fprintf('Data label: %s\n',data_label);
fprintf('Filter order: %d\n',order);