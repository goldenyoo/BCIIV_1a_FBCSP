function [result] = myClassifier(evaluation_trial,training_data,training_label)
    a= 100;
    b=100;
    
    p_x_w_1 = 1;
    p_x_w_2 = 1;
    for j = 1:length(evaluation_trial)
        xj = evaluation_trial(j);
        f_j = training_data(:,j);
        tmp1 = xj - f_j(find(training_label == 1));
        tmp2 = xj - f_j(find(training_label == -1));
        
        ha = ((4/(3*a))^0.2)*(var(tmp1))^0.5;
        hb = ((4/(3*b))^0.2)*(var(tmp2))^0.5;
        
        pi_1 = (1/sqrt(2*pi))*exp((-1/(2*ha*ha))*tmp1.*tmp1);
        pi_2 = (1/sqrt(2*pi))*exp((-1/(2*hb*hb))*tmp2.*tmp2);
        
        p_hat_xj_w_1 = (1/a)*sum(pi_1);
        p_hat_xj_w_2 = (1/b)*sum(pi_2);
        
        p_x_w_1 = p_x_w_1*p_hat_xj_w_1;
        p_x_w_2 = p_x_w_2*p_hat_xj_w_2;
    end
    
    MAP_1 = 0.5*p_x_w_1/(0.5*p_x_w_1+0.5*p_x_w_2);
    MAP_2 = 0.5*p_x_w_2/(0.5*p_x_w_1+0.5*p_x_w_2);
 
    if MAP_1 > MAP_2
        result = 1;  
    else
        result = -1;   % r
    end
end
