% ----------------------------------------------------------------------- %
%    File_name: Eval.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_02_10                           
%                                                            
 % ----------------------------------------------------------------------- %
function [Sc, Mse] = Eval(answer,interest_freq_band,interest_P, training_data,training_label,ref)
data_label = string(answer(1,1));   
m = double(string(answer(2,1))); % feature vector will have length (2m)
referencing = double(string(answer(3,1))); % Non(0), CAR(1), LAP(2)
order = double(string(answer(4,1))); % Filter order


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

FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_eval_ds1',data_label,'.mat');
chunk = 350;
fs = 100;

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
        cnt_y = cnt(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)        
        Means = (1/size(cnt_y,1))*sum(cnt_y);
        for i = 1 : size(cnt_y,1)
            cnt_y(i,:) = cnt_y(i,:) - Means; % CAR
        end
     % LAP   
    elseif referencing == 2 
        cnt_n = myLAP(cnt,nfo); % Laplacian
        cnt_y = cnt_n(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
    end
else
        %%% Calculate differential voltage
    for i = 1 : size(cnt,1)
        cnt(i,:) = cnt(i,:) - cnt(ref,:);
    end
    
    cnt_y = cnt(3:55,:); % Exclude electrode (AF3, AF4, O1, O2, PO1, PO2)
end

clear cnt 

%% 
for fb = 1:size(interest_freq_band,1)
        low_f = interest_freq_band(fb,1);
        high_f = interest_freq_band(fb,2);
        %BPF Design
        bpFilt = designfilt('bandpassfir','FilterOrder',order, ...
            'CutoffFrequency1',low_f,'CutoffFrequency2',high_f, ...
            'SampleRate',fs);
        
%          bpFilt = designfilt('bandpassiir','FilterOrder',order, ...
%     'StopbandFrequency1',low_f,'StopbandFrequency2',high_f, ...
%     'StopbandAttenuation',40,'SampleRate',fs);

        
        % Apply BPF
        for i = 1:size(cnt_y,1)
            cnt_x(i,:) = filtfilt(bpFilt, cnt_y(i,:));            
        end
        filtered{fb} = cnt_x;
end

%% 

score = [];
predictions = [];
checks = [];

% For class 2
nbytes = 0;
for j = 1 : length(B)
    tmp1 = round(B(j,1)/10);
    tmp2 = round(B(j,2)/10);
    
%     fprintf(repmat('\b',1,nbytes));
%     nbytes = fprintf('class 2 --> %d / %d ',j,length(B));
    
    for fb = 1:size(interest_freq_band,1)
        cnt_c = filtered{fb};
                
        E = cnt_c(:, tmp1:tmp2);
        Z = interest_P{fb}'*E;
        
        % Feature vector
        tmp_ind = size(Z,1);
        Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];

        var_vector = diag(Z_reduce*Z_reduce')/trace(Z_reduce*Z_reduce');
        fp = log(var_vector);
        
        evaluation_trial(1,1+2*m*(fb-1):2*m*fb) = fp;    
    end 
       
    % Run classifier
    
    [prediction] = myClassifier(evaluation_trial,training_data,training_label);
    if prediction == -1
        score = [score 1];
        
    else
        score = [score 0];
        
    end
    predictions = [predictions prediction];
       
end


nbytes = 0;
% For class 1
for j = 1 : length(D)
    tmp1 = round(D(j,1)/10);
    tmp2 = round(D(j,2)/10);
        
    for fb = 1:size(interest_freq_band,1)
        cnt_c = filtered{fb};
                
        E = cnt_c(:, tmp1:tmp2);
        Z = interest_P{fb}'*E;
        
        % Feature vector
        tmp_ind = size(Z,1);
        Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];

        var_vector = diag(Z_reduce*Z_reduce')/trace(Z_reduce*Z_reduce');        
        fp = log(var_vector);
        
        evaluation_trial(1,1+2*m*(fb-1):2*m*fb) = fp;    
    end   
    
      
    % Run classifier
    [prediction] = myClassifier(evaluation_trial,training_data,training_label);
    if prediction == 1
        score = [score 1];
        
    else
        score = [score 0];
        
    end
    predictions = [predictions prediction];
       
end


% Caculation score
Sc = 100*sum(score)/length(score);
Mse =  length(find(score == 0))*4/length(score);
end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
