% ----------------------------------------------------------------------- %
%    File_name: Eval.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_27                           
%                                                            
 % ----------------------------------------------------------------------- %
function output = Eval(answer,Mr,Ml,Qr,Ql,P,ref)
data_label = string(answer(1,1));   
m = double(string(answer(2,1))); % feature vector will have length (2m)
low_f = double(string(answer(3,1))); % Low cutoff freq
high_f = double(string(answer(4,1))); % High cutoff freq
sampling_rate = double(string(answer(5,1)));
referencing = double(string(answer(6,1))); % Non(0), CAR(1), LAP(2)
order = double(string(answer(7,1))); % Filter order


%% Call true label
FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\true_labels\BCICIV_eval_ds1',data_label,'_1000Hz_true_y.mat');
load(FILENAME);

A=[];
B=[];
C=[];
D=[];
%NaN, -1, 0, 1
tmp = 1;
for i=1:length(true_y)-1
    if isnan(true_y(i,1))
        if isnan(true_y(i+1,1))
            continue
        else
            A = [A; [tmp i]];
            tmp = i+1;
        end
    elseif true_y(i,1) ~= true_y(i+1,1)
        if true_y(i,1)== -1
            B = [B; [tmp i]];
        elseif true_y(i,1)== 0
            C = [C; [tmp i]];
        elseif true_y(i,1)== 1
            D = [D; [tmp i]];
        end
        tmp = i+1;
    end
end

clear true_y
%% 
% Load file
if sampling_rate == 0
    FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_eval_ds1',data_label,'.mat');
    chunk = 350;
    fs = 100;
else
    FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1eval_1000Hz_mat\BCICIV_eval_ds1',data_label,'_1000Hz.mat');
    chunk = 3500;
    fs = 1000;
end
load(FILENAME);

% Data rescale
cnt= 0.1*double(cnt);
cnt = cnt';
%% Preprocessing
if referencing ~= 0
    %%% Calculate differential voltage
    for i = 1 : size(cnt,1)
        cnt(i,:) = cnt(i,:) - cnt(ref,:);
    end

    % common average
    if referencing == 1         
        cnt_c = cnt(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)        
        Means = (1/size(cnt_c,1))*sum(cnt_c);
        for i = 1 : size(cnt_c,1)
            cnt_c(i,:) = cnt_c(i,:) - Means; % CAR
        end
     % LAP   
    elseif referencing == 2 
        cnt_n = myLAP(cnt,nfo); % Laplacian
        cnt_c = cnt_n(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
    end
else
        %%% Calculate differential voltage
    for i = 1 : size(cnt,1)
        cnt(i,:) = cnt(i,:) - cnt(ref,:);
    end
    
    cnt_c = cnt(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
end

clear cnt cnt_n

%% 
%BPF Design
bpFilt = designfilt('bandpassfir','FilterOrder',order, ...
         'CutoffFrequency1',low_f,'CutoffFrequency2',high_f, ...
         'SampleRate',fs);

% Apply BPF
for i = 1:size(cnt_c,1)
    cnt_c(i,:) = filtfilt(bpFilt, cnt_c(i,:));
end

%% 
% f1 = figure;
% f2 = figure;
score = [];
predictions = [];
checks = [];

% For class 1
for j = 1 : length(B)
    if sampling_rate == 0
        tmp1 = round(B(j,1)/10);
        tmp2 = round(B(j,2)/10);
    else
        tmp1 = B(j,1);
        tmp2 = B(j,2);
    end
    E = cnt_c(:, tmp1:tmp2);
    Z = P'*E;
    
     % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
   
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
        
    fp = log(var_vector);
    fp = fp';
    
    % Graphical represent
%      figure(f1)
%      scatter3(Z(1,:), Z(size(cnt_c,1),:),Z(2,:),'b'); hold on;
%      figure(f2)
%      scatter3(fp(1),fp(2),fp(6),'b'); hold on;
       
    % Run classifier
    [check, prediction] = myClassifier(fp,Mr,Ml,Qr,Ql);
    if prediction == -1
        score = [score 1];
        
    else
        score = [score 0];
        
    end
    predictions = [predictions prediction];
    checks = [checks check];
    
end

% For class 2
for j = 1 : length(D)
    if sampling_rate == 0
        tmp1 = round(D(j,1)/10);
        tmp2 = round(D(j,2)/10);
    else
        tmp1 = D(j,1);
        tmp2 = D(j,2);
    end
    E = cnt_c(:, tmp1:tmp2);
    Z = P'*E;
    
    % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
   
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
        
    fp = log(var_vector);
    fp = fp';
    
    % Graphical represent
%      figure(f1)
%      scatter3(Z(1,:), Z(size(cnt_c,1),:),Z(2,:),'r'); hold on;
%      figure(f2)
%      scatter3(fp(1),fp(2),fp(6),'r'); hold on;
       
    % Run classifier
    [check, prediction] = myClassifier(fp,Mr,Ml,Qr,Ql);
    if prediction == 1
        score = [score 1];
        
    else
        score = [score 0];
        
    end
    predictions = [predictions prediction];
    checks = [checks check];
    
end

% Caculation score
output = 100*sum(score)/length(score);

end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
