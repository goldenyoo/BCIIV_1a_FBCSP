% ----------------------------------------------------------------------- %
%    File_name: FBCSP.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_02_18                            
%                                                            
 % ----------------------------------------------------------------------- %

%% Get input parameter from user
close all
clear all

% Ask user for input parameters
prompt = {'Data label: ', 'Feature vector length: '};
dlgtitle = 'Input';
dims = [1 50];
definput = {'a', '2'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
% Error detection
if isempty(answer), error("Not enough input parameters."); end

%% Conditions
% Rereferencing method 
ref_method = [0 1 2]; % Non(0), CAR(1), LAP(2)

% Filter order
filt_ord = [14]; %%%%%%%%%%%%%%%%%%%%%% Change

% Reference electrode number
ref = 29;        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change

%% CSP 
for i = 1:length(ref_method)
    for j = 1:length(filt_ord)
        answer(3,1) = {ref_method(i)};
        answer(4,1) = {filt_ord(j)};
        fprintf('Data_Label: %s\n',string(answer(1,1)));
        fprintf('Re-referencing: %d',ref_method(i));
        fprintf(' filter_order: %d\n',filt_ord(j));
        [interest_freq_band,interest_P, training_data,training_label] = Calib(answer,ref);
        
        for k = 1:size(interest_freq_band,1)
            fprintf('Filter bank: %d %d\n',interest_freq_band(k,1),interest_freq_band(k,2));
        end
        
        [score, mse] = Eval(answer,interest_freq_band,interest_P, training_data,training_label,ref);
        
        fprintf(' ----> score: %f\n\n',score);
%         fprintf(' ----> mse: %f\n\n',mse);
         output(i,j) = score;
    end
    fprintf('\n');
end
disp(output);
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
