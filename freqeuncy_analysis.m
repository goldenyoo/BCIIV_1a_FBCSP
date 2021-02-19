%% Call raw data
close all
clear all

% Ask user for input parameters
prompt = {'Data label: '};
dlgtitle = 'Input';
dims = [1 50];
definput = {'a'};
answer = inputdlg(prompt,dlgtitle,dims,definput);


% Error detection
if isempty(answer), error("Not enough input parameters."); end

% Input parameters
data_label = string(answer(1,1));   % Calib_ds1 + "data_label"



ref=29;


% Load file
FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_calib_ds1',data_label,'.mat');
% FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_eval_ds1',data_label,'.mat');
load(FILENAME);

% Data rescale
cnt= 0.1*double(cnt);
cnt = cnt';

%% Preprocessing

%%% Calculate differential voltage
for i = 1 : size(cnt,1)
    cnt_k(i,:) = cnt(i,:) - cnt(ref,:);
end

% common average

cnt_y = cnt_k(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
Means = (1/size(cnt_y,1))*sum(cnt_y);
for i = 1 : size(cnt_y,1)
    cnt_y(i,:) = cnt_y(i,:) - Means; % CAR
end
cnt_yn(3:55,:) = cnt_y;
cnt_yn(56:59,:) = 0;


cnt_n = myLAP(cnt,nfo); % Laplacian

   bpFilt = designfilt('bandpassiir','SampleRate',100,'PassbandFrequency1',8, ...
    'PassbandFrequency2',12,'StopbandFrequency1',8-2,'StopbandFrequency2',12+2, ...
    'StopbandAttenuation1',60,'StopbandAttenuation2',60, 'PassbandRipple',1,'DesignMethod','cheby2');

% Apply BPF
for i = 1:size(cnt_n,1)
    cnt_x_CAR(i,:) = filtfilt(bpFilt, cnt_yn(i,:));
end

% Apply BPF
for i = 1:size(cnt_n,1)
    cnt_x_NON(i,:) = filtfilt(bpFilt, cnt_k(i,:));
end

% Apply BPF
for i = 1:size(cnt_n,1)
    cnt_x_LAP(i,:) = filtfilt(bpFilt, cnt_n(i,:));
end

%% 
%%%%%% cnt, cnt_k, cnt_yn, cnt_n

% for i = [27 31] 
for i = 3:55
%     time = 0:1/100:(190594-1)/100; 
    
%     figure
%     subplot(4,1,1)
%     sinad(cnt_k(i,:),100);
%     title("Differential")
%     legend('off')
%     
%     subplot(4,1,2)
%     sinad(cnt_x_10(i,:),100);
%     title("Order: 10")
%     legend('off')
%     
%     subplot(4,1,3)
%     sinad(cnt_x_15(i,:),100);
%     title("Order: 15")
%     legend('off')
%     
%     subplot(4,1,4)
%     sinad(cnt_x_30(i,:),100);
%     title("Order: 30")
%     legend('off')

    figure
    subplot(4,1,1)
    sinad(cnt(i,:),100);
    title("Raw")
    legend('off')
    
    subplot(4,1,2)
    sinad(cnt_x_NON(i,:),100);
    title("NON")
    legend('off')
    
    subplot(4,1,3)
    sinad(cnt_x_CAR(i,:),100);
    title("CAR")
    legend('off')
    
    subplot(4,1,4)
    sinad(cnt_x_LAP(i,:),100);
    title("LAP")  
    legend('off')
    
    sgtitle(nfo.clab{i});
    
end