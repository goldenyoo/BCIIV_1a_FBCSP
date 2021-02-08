% ----------------------------------------------------------------------- %
%    File_name: Calib.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_27                           
%                                                            
 % ----------------------------------------------------------------------- %
function [Mr,Ml,Qr,Ql,P] = Calib(answer,ref)

% Input parameters
data_label = string(answer(1,1));   
m = double(string(answer(2,1))); % feature vector will have length (2m)
low_f = double(string(answer(3,1))); % Low cutoff freq
high_f = double(string(answer(4,1))); % High cutoff freq
sampling_rate = double(string(answer(5,1)));
referencing = double(string(answer(6,1))); % Non(0), CAR(1), LAP(2)
order = double(string(answer(7,1))); % Filter order

% Load file
if sampling_rate == 0
    FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1_mat\BCICIV_calib_ds1',data_label,'.mat');
    chunk = 350;
    fs = 100;
else
    FILENAME = strcat('C:\Users\유승재\Desktop\Motor Imagery EEG data\BCICIV_1calib_1000Hz_mat\BCICIV_calib_ds1',data_label,'_1000Hz.mat');
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

%% Calculate spatial filter
a = 1; b = 1;
C_r = zeros(size(cnt_c,1)); C_l = zeros(size(cnt_c,1));

% Training only for training data set
for i = 1:length(mrk.pos)
    
    % One trial data
    E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+chunk);
    
    %Centering
    %     Means = mean(E,2);
    %     E = E - diag(Means)*ones(size(E,1),size(E,2));
        
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
C_r = C_r/(a-1);
C_l = C_l/(b-1);

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
[U, phsi] = eig(Sa);

% sort
[d, ind] = sort(abs(diag(phsi)),'descend');
phsi_new = diag(d);
U_new = U(:,ind);

% Total Projection matrix,   Z = P'*X
P = (U_new'*W)';

%% Calculate feature vector

fp_r = [];
fp_l = [];

for i = 1:length(mrk.pos)
    
    % One trial data
    E = cnt_c(:,mrk.pos(1,i):mrk.pos(1,i)+chunk);
    
    %     % Centering
    %     Means = mean(E,2);
    %     E = E - diag(Means)*ones(size(E,1),size(E,2));    
    
    % Project data using calculated spatial filter
    Z = P'*E;
    
    % Feature vector
    tmp_ind = size(Z,1);
    Z_reduce = [Z(1:m,:); Z(tmp_ind-(m-1):tmp_ind,:)];
    var_vector = var(Z_reduce,0,2)';
    var_vector = (1/sum(var_vector))*var_vector;
    
    fp = log(var_vector);
    fp = fp';
    
    if mrk.y(1,i) == 1
        fp_r = [fp_r fp];
    else
        fp_l = [fp_l fp];
    end
end

Mr = mean(fp_r,2);
Ml = mean(fp_l,2);

Qr = zeros(2*m);
for i = 1:length(fp_r)
    tmp = (fp_r(:,i) - Mr)*(fp_r(:,i) - Mr)';
    Qr = Qr + tmp;
end
Qr = (1/(length(fp_r)-1))*Qr;

Ql = zeros(2*m);
for i = 1:length(fp_l)
    tmp = (fp_l(:,i) - Ml)*(fp_l(:,i) - Ml)';
    Ql = Ql + tmp;
end
Ql = (1/(length(fp_l)-1))*Ql;

end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
