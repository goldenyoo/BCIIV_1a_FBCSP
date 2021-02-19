% ----------------------------------------------------------------------- %
%    File_name: Calib.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_02_18                           
%                                                            
 % ----------------------------------------------------------------------- %
function [interest_freq_band,interest_P, training_data,training_label] = Calib(answer,ref)

% Input parameters
data_label = string(answer(1,1));
m = double(string(answer(2,1))); % feature vector will have length (2m)
referencing = double(string(answer(3,1))); % Non(0), CAR(1), LAP(2)
order = double(string(answer(4,1))); % Filter order

% Load file

FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_calib_ds1',data_label,'.mat');
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
% FB = [[4 8];[10 14]; [12 16]; [16 20]; [20 24]; [24 28]; [28 32]; [32 36]; [36 40]];
% FB = [[4 8];[8 12]; [12 16]; [16 20]];
FB = [[8 10]; [10 12]; [12 14]; [14 16]];

feature_V = [];

for fb = 1:length(FB)
    low_f = FB(fb,1);
    high_f = FB(fb,2);
    
    %BPF Design
%     bpFilt = designfilt('bandpassfir','FilterOrder',order, ...
%         'CutoffFrequency1',low_f,'CutoffFrequency2',high_f, ...
%         'SampleRate',fs);

%      bpFilt = designfilt('bandpassiir','FilterOrder',order, ...
%     'StopbandFrequency1',low_f,'StopbandFrequency2',high_f, ...
%     'StopbandAttenuation',40,'SampleRate',fs);
    
    bpFilt = designfilt('bandpassiir','SampleRate',fs,'PassbandFrequency1',low_f, ...
    'PassbandFrequency2',high_f,'StopbandFrequency1',low_f-2,'StopbandFrequency2',high_f+2, ...
    'StopbandAttenuation1',40,'StopbandAttenuation2',40, 'PassbandRipple',1,'DesignMethod','cheby2');
    
    % Apply BPF
    for i = 1:size(cnt_y,1)
        cnt_x(i,:) = filtfilt(bpFilt, cnt_y(i,:));
    end
    filtered{fb} = cnt_x;
end


%% 

for fb = 1:length(FB)
  
    cnt_c = filtered{fb};
    
    % Calculate spatial filter
    a = 0; b = 0;
    C_r = zeros(size(cnt_c,1)); C_l = zeros(size(cnt_c,1));
   
    % Training only for training data set
    for i = 1:length(mrk.pos)
        
        % One trial data
        E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+chunk);        
        
        % Covariance 연산
        C = E*E'/ trace( E*E');
        
        % According to its class, divide calculated covariance
        if mrk.y(1,i) == 1
            C_r = C_r+C;
            a = a+1;            
        else
            C_l = C_l+C;
            b = b+1;
        end
        
    end
    
    % Average covariance of each class
    C_r = C_r/(a);
    C_l = C_l/(b);
    
    % composite covariace
    C_c = C_r + C_l;
    
    % EVD for composite covariance
    [V, D] = eig(C_c);
    
    % sort eigen vector with descend manner
    [d, ind] = sort(abs(diag(D)),'descend');
    D_new = diag(d);
    V_new = V(:,ind);
    
    % whitening transformation
    whiten_tf = V_new*D_new^(-0.5);
    W = whiten_tf';
    
    % Apply whitening to each averaged covariances
    Sa = W*C_r*W';
    Sb = W*C_l*W';
    
    % EVD for transformed covariance
    [U, phsi] = eig(Sa,Sa+Sb);
    
    % sort
    [d, ind] = sort(abs(diag(phsi)),'descend');
    phsi_new = diag(d);
    U_new = U(:,ind);
    
    % Total Projection matrix,   Z = P'*X
    P = (U_new'*W)';
    P_FB{fb} = P;    
    
    % Calculate feature vector
    for i = 1:length(mrk.pos)
        
        % One trial data
        E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+chunk);
        
        % Project data using calculated spatial filter
        Z = P'*E;
        
        % Feature vector
        tmp_ind = size(Z,1);
        Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
%         var_vector = var(Z_reduce,0,2)';
%         var_vector = (1/sum(var_vector))*var_vector;
        
        var_vector = diag(Z_reduce*Z_reduce')/trace(Z_reduce*Z_reduce');
        fp = log(var_vector);
        
        feature_V(i,1+2*m*(fb-1):2*m*fb) = fp;        
    end
end



% Caculate Mutual information
for j = 1:size(feature_V,2)
    f_j = feature_V(:,j);
    H_w = -(0.5*log2(0.5)+0.5*log2(0.5));
    
    p_fj_w_1 = 1;  p_fj_w_2 = 1;
    for i = 1:size(feature_V,1)
        
        tmp1 = f_j(i) - f_j(find(mrk.y(1,:) == 1));
        tmp2 = f_j(i) - f_j(find(mrk.y(1,:) == -1));
        
        ha = ((4/(3*a))^0.2)*(var(tmp1))^0.5;
        hb = ((4/(3*b))^0.2)*(var(tmp2))^0.5;
        
        pi_1 = (1/sqrt(2*pi))*exp((-1/(2*ha*ha))*tmp1.*tmp1);
        pi_2 = (1/sqrt(2*pi))*exp((-1/(2*hb*hb))*tmp2.*tmp2);
        
        p_hat_fji_w_1 = (1/a)*sum(pi_1);
        p_hat_fji_w_2 = (1/b)*sum(pi_2);
        
        p_fj_w_1 = p_fj_w_1*p_hat_fji_w_1;
        p_fj_w_2 = p_fj_w_2*p_hat_fji_w_2;
    end
    
    p_w_1_fj = p_fj_w_1*0.5/(p_fj_w_1*0.5 + p_fj_w_2*0.5);
    p_w_2_fj = p_fj_w_2*0.5/(p_fj_w_1*0.5 + p_fj_w_2*0.5);    
    
    H_w_fj = -(p_w_1_fj*log2(p_w_1_fj) + p_w_2_fj*log2(p_w_2_fj));
    
    if isnan(H_w_fj)
        H_w_fj = 0;
    end
    
    H_w_f(1,j) = H_w_fj;
    MI_j = H_w - H_w_fj;
    
    MI(1,j) = MI_j;
end


[d, ind_tmp] = sort(H_w_f);
ind = ind_tmp(1:2); %%%%%%%%%%%%%%%% Change number of k %%%%%%%%%%%%%

peter = [];
for k = ind
    if fix((k-1)/(2*m))==0
        peter = [peter 1];
    elseif fix((k-1)/(2*m))==1
        peter = [peter 2];
    elseif fix((k-1)/(2*m))==2
        peter = [peter 3];
    elseif fix((k-1)/(2*m))==3
        peter = [peter 4];
    elseif fix((k-1)/(2*m))==4
        peter = [peter 5];
    elseif fix((k-1)/(2*m))==5
        peter = [peter 6];
    elseif fix((k-1)/(2*m))==6
        peter = [peter 7];
    elseif fix((k-1)/(2*m))==7
        peter = [peter 8];
    elseif fix((k-1)/(2*m))==8
        peter = [peter 9];
    end        
end

% press = unique([peter]);
% press = unique([peter 2]); %%%%%%%%%%%%%% If you want 8~12Hz always included
press = [1 2];

interest_freq_band = FB(press,:);
interest_P = {};


for pt=1:length(press)
    interest_P{pt} = P_FB{press(pt)};
end
training_data=[];
for pt = press
    training_data = [training_data feature_V(:,1+2*m*(pt-1):2*m*(pt))];
    
end
training_label = mrk.y(1,:);
end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
