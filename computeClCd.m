function [cl,cd]=computeClCd(alpha,M)
    %Uses data extracted from Xflr5
    %At Re=2.4e6 and M=0.5;
    
    load Data_xflr_file;
    
    cl=interp1q(aClaCd(:,1),aClaCd(:,2),alpha);
    cd=interp1q(aClaCd(:,3),aClaCd(:,4),alpha);
    
    cl=cl/(sqrt(1-M^2));
    cd=cd/(sqrt(1-M^2));

end