function [dcl,dcd]=compute_derClCd(alpha,M)

    %Uses data extracted from Xflr5
    %At Re=2.4e6 and M=0;

    load derivades.mat;
   
    dcl=interp1q(adCladCd(:,1),adCladCd(:,2),alpha);
    dcd=interp1q(adCladCd(:,3),adCladCd(:,4),alpha);
    
    dcl=dcl/(sqrt(1-M^2));
    dcd=dcd/(sqrt(1-M^2));

end