% ----------------------------------------------------------------------- %
%    File_name: CSP.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_27                            
%                                                            
 % ----------------------------------------------------------------------- %

%% Get input parameter from user
close all
clear all

% Ask user for input parameters
prompt = {'Data label: ', 'Feature vector length: '};
dlgtitle = 'Input';
dims = [1 50];
definput = {'a', '3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
% Error detection
if isempty(answer), error("Not enough input parameters."); end

%% Conditions
% Rereferencing method 
ref_method = [0]; % Non(0), CAR(1), LAP(2)

% Filter order
filt_ord = [25 ];

% Reference electrode number
ref = 29;        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change

%% CSP 
for i = 1:length(ref_method)
    for j = 1:length(filt_ord)
        answer(3,1) = {ref_method(i)};
        answer(4,1) = {filt_ord(j)};
        fprintf('Re-referencing: %d',ref_method(i));
        fprintf(' filter_order: %d',filt_ord(j));
        [Mr,Ml,Qr,Ql,P] = Calib(answer,ref);
        tmp = Eval(answer,Mr,Ml,Qr,Ql,P,ref);
        output(i,j) = tmp;
        fprintf(' ----> score: %f\n',tmp);
    end
    fprintf('\n');
end

%% Output present

fprintf('Data_Label: %s\n',string(answer(1,1)));
fprintf('BPF cutoff freq: %s ~ %s\n',string(answer(3,1)),string(answer(4,1)));
disp(output);
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %
