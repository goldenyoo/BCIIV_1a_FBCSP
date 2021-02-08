function [check, result] = myClassifier(fp,Mr,Ml,Qr,Ql)
    
      
    check = (0.5*(fp - Mr)'*pinv(Qr)*(fp - Mr) + 0.5*log(det(Qr))) - (0.5*(fp - Ml)'*pinv(Ql)*(fp - Ml) + 0.5*log(det(Ql))); 
 
    if check > 0
        result = -1;  % l
    else
        result = 1;   % r
    end
end
