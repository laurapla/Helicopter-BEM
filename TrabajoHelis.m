clear all;
clc;

m=4500;                
Cl = @(alpha)  1E-06.*alpha.^4-0.0001.*alpha.^3-0.0008.*alpha.^2+0.1197.*alpha+0.259;  
Cd = @(alpha)  4E-07.*alpha.^4 - 3E-06.*alpha.^3 - 1E-05.*alpha.^2 + 0.001.*alpha + 0.0056;  
alphaOPT=5;   %XFLR5        

g=9.81;        
W=m*g;         
DL=35;         
To=W;          
A=To/(g*DL);    
R=(A/pi)^0.5;  
h=1500;         
Mtip=0.5;      
nb=3;           %Initial
P0=101325; 
P=84556;        %ISA
rho0=1.225; 
rho=1.05807;    %ISA
T0=288.15; 
T=278.400;      %ISA
a0=340.29;
a=334.487;      %ISA
Omega=a*Mtip/R;   %Omega=a*Mtip/R

%Momenthum theory

Vc=0;                      
Vio=(To/(2*rho*A))^0.5;    
lambda=Vio/(Omega*R);  
alpha=alphaOPT;            
theta = @(alpha,lambda,r) alpha*(pi/180)+atan(lambda/r);
sigma = @(alpha,lambda,r) (8*r*lambda^2)/((r^2+lambda^2)*(Cl(alpha)*cos(theta(alpha,lambda,r)-alpha*(pi/180))-Cd(alpha)*sin(theta(alpha,lambda,r)-alpha*(pi/180))));
for i=1:1:21
    r(i)=0.01*(i+60)-0.01;     
    t_07(i)=theta(alpha,lambda,r(i));
    s_07(i)=sigma(alpha,lambda,r(i));
   
end
p = polyfit(r,s_07,1);         
sigmalin = @(r) p(1)*r+p(2); 
q = polyfit(r,t_07,1);         
thetalin = @(r) q(1)*r+q(2);  
for i=1:1:101
    r(i)=0.01*i-0.01;                  
    s(i)=sigma(alpha,lambda,r(i));  
    s_lin(i)=sigmalin(r(i));          
    t(i)=theta(alpha,lambda,r(i));                   
    tlin(i)=thetalin(r(i));         
    lambda_i_MTH(i)=lambda;        
end 

Croot=(sigmalin(0)*pi*R)/nb;  
if Croot < 0.5                
    nb;                      
else                          
    while Croot > 0.5         
    nb=nb+1;    
    Croot=(sigmalin(0)*pi*R)/nb;
    end
    
    if nb<7
        nb; 
    else   
        nb=3;
        while Croot < 0.75    
            nb=nb+1;    
            Croot=(sigmalin(0)*pi*R)/nb;
        end
    nb; 
    end        
end


for i=1:1:101 
    c(i)=(sigmalin(r(i))*pi*R)/nb;
end
 
%BEM

for i=1:1:6
Vc(i)=(i-1)*2.5;  
lc1(i)=Vc(i)/(Omega*R); 
syms xvi
eqn = To == 2*rho*A*xvi*(xvi+Vc(i));  
Vinduced=solve(eqn,xvi);
Vi(i)=abs(Vinduced(1));
end


for k=1:1:6
j=0;
T=0;
t0=0.1;
    while T < To
    t0=t0+j*0.0001;      
        for i=1:1:101      
            tlin(i)=q(1)*r(i)+t0;  
        end  
            for i=1:1:101    
            lc=lc1(k);
            radius=r(i);
            sig=sigmalin(radius);
            tor=tlin(i);
            BEM=@(li) (8.*((li+lc).*li).*radius)-(sig.*(radius^2+(li+lc).^2).*(Cl((180/pi).*(tor-atan((li+lc)./radius))).*cos(atan((li+lc)./radius))-Cd((180/pi).*(tor-atan((li+lc)./radius))).*sin(atan((li+lc)./radius))));          
            l_i(i) = fzero(BEM,0.025);         
           
            end
    
        
        for i=1:1:101
            y(i)=r(i)*R;                               
            dy=R/100;                                  
            U=sqrt((Omega*y(i))^2+((lc1(k)+l_i(i))*Omega*R)^2);  
            phi=atan((l_i(i)+lc1(k))/r(i));      
            alpha=(180/pi)*(tlin(i)-phi);            
            dL(i)=0.5*rho*(U^2)*c(i)*dy*Cl(alpha);     
            dD(i)=0.5*rho*(U^2)*c(i)*dy*Cd(alpha);
            dFz(i)=dL(i)*cos(phi)-dD(i)*sin(phi);
            OmegaydLx(i)=Omega*y(i)*dL(i)*sin(phi);   
            OmegaydDx(i)=Omega*y(i)*dD(i)*cos(phi);   
        end
        T=sum(dFz)*nb;        
    j=j+1;
    end
    thrust(k)=T;
    pitch(k)=t0;

%POWER

Pi_1(k)=2*rho*A*Vi(k)*(Vi(k)+Vc(k))^2;  % Induced power using MTH
Po1v=0.5*rho*(Omega*R)^3*pi*R^2*s.*(r.^3);
Po_1=trapz(Po1v);
Pi_2(k)=nb*sum(OmegaydLx);        % Induced power using BEM
Po_2(k)=nb*sum(OmegaydDx);        % Parasite power using BEM
PT_2(k)=Pi_2(k)+Po_2(k);      % Total power using BEM

% lift = polyfit(y,dL,3);         
% dLift = @(y) lift(1)*y.^3+lift(2)*y.^2+lift(3)*y+lift(4);
% dr = polyfit(r,dD,3);         
% dDrag = @(r) dr(1)*r.^3+dr(2)*r.^2+dr(3)*r+dr(4);
% for i=1:1:101 
%     dL2(k,i)=dL(i);
%     dD2(k,i)=dD(i);
%     dLpolinomial(k,i)=dLift(y(i));    
%     dDpolinomial(k,i)=dDrag(r(i)); 
%     lambdai2(k,i)=l_i(i);
%     tlin2(k,i)=tlin(i);
% end

end

% for i=1:1:101 %the subindex indicates the ascensional velocity
% 
%  dL0(i)=dL2(1,i);    
%  dD0(i)=dD2(1,i);
%  dL25(i)=dL2(2,i);    
%  dD25(i)=dD2(2,i);
%  dL50(i)=dL2(3,i);    
%  dD50(i)=dD2(3,i);
%  dL75(i)=dL2(4,i);    
%  dD75(i)=dD2(4,i);
%  dL100(i)=dL2(5,i);    
%  dD100(i)=dD2(5,i);
%  dL125(i)=dL2(6,i);    
%  dD125(i)=dD2(6,i);
%  
%  dLpol0(i)=dLpolinomial(1,i);    
%  dDpol0(i)=dDpolinomial(1,i);
%  dLpol25(i)=dLpolinomial(2,i);    
%  dDpol25(i)=dDpolinomial(2,i);
%  dLpol50(i)=dLpolinomial(3,i);    
%  dDpol50(i)=dDpolinomial(3,i);
%  dLpol75(i)=dLpolinomial(4,i);    
%  dDpol75(i)=dDpolinomial(4,i);
%  dLpol100(i)=dLpolinomial(5,i);    
%  dDpol100(i)=dDpolinomial(5,i);
%  dLpol125(i)=dLpolinomial(6,i);    
%  dDpol125(i)=dDpolinomial(6,i);
%  
%  lambdai1(i)=lambdai2(1,i);
%  tlin1(i)=tlin2(1,i);
% end
% 
% T_BEM_0=thrust(1)             
% Pitch_BEM_0=pitch(1)*(180/pi) 
% T_BEM_2dot5=thrust(2)             
% Pitch_BEM_2dot5=pitch(2)*(180/pi)
% T_BEM_5=thrust(3)             
% Pitch_BEM_5=pitch(3)*(180/pi)
% T_BEM_7dot5=thrust(4)             
% Pitch_BEM_7dot5=pitch(4)*(180/pi)
% T_BEM_10=thrust(5)             
% Pitch_BEM_10=pitch(5)*(180/pi)
% T_BEM_12dot5=thrust(6)             
% Pitch_BEM_12dot5=pitch(6)*(180/pi)


% figure    
%     subplot(2,2,1);
%     plot(r,s,'r',r,s_lin,'b');
%     xlabel('r');ylabel('\sigma');
%     title('Solidity');
%     legend('Ideal','Linear','Location','best');grid;
%     
%     subplot(2,2,3);
%     plot(r,t,'r',r,tlin1,'b');
%     xlabel('r');ylabel('\theta');
%     title('\theta');
%     legend('Ideal','Linear','Location','best');grid;
%     
%     subplot(2,2,2); 
%     plot(r,c,'b');
%     xlabel('r');ylabel('c');
%     title('Blade chord distribution');grid;
%         
%     subplot(2,2,4); 
%     plot(r,lambdai1,'r',r,lambda_i_MTH,'b');
%     xlabel('r');ylabel('\lambda_i');
%     title('Induced velocities distribution');
%     legend('BEM','MTH','Location','best');grid; 
%     
%     
% figure
%    
%     subplot(2,1,1);
%     plot(Vc,Pi_1,'r',Vc,Pi_2,'b');
%     xlabel('Vc');ylabel('W');
%     title('Velocity vs Induced Power');
%     legend('MTH','BEM','Location','best');grid; 
%     
%    
%     subplot(2,1,2);
%     plot(Vc,Pi_2,'r',Vc,Po_2,'b',Vc,PT_2,'g');
%     xlabel('Vc');ylabel('W');
%     title('Velocity vs Power');
%     legend('Induced power','Parasite power','Total power','Location','best');grid;
%     
% figure
%     
%     subplot(3,2,1);
%     plot(y,dL0,'r-o',y,dLpol0,'b');
%     xlabel('y');ylabel('N/m');title('Lift Distribution 0m/s');
%    
%     
%     subplot(3,2,3); 
%     plot(y,dL25,'r-o',y,dLpol25,'b');
%     xlabel('y');ylabel('N/m');title('Lift Distribution 2,5m/s');
%  
%     
%     subplot(3,2,5);
%     plot(y,dL50,'r-o',y,dLpol50,'b');
%     xlabel('y');ylabel('LN/m');title('Lift Distribution 5m/s');
%    
%     
%     subplot(3,2,2);
%     plot(y,dL75,'r-o',y,dLpol75,'b');
%     xlabel('y');ylabel('N/m');title('BLift Distribution 7,5m/s');
%    
%     
%     subplot(3,2,4);
%     plot(y,dL100,'r-o',y,dLpol100,'b');
%     xlabel('y');ylabel('N/m');title('Lift Distribution 10m/s');
%     
%     
%     subplot(3,2,6);
%     plot(y,dL125,'r-o',y,dLpol125,'b');
%     xlabel('y');ylabel('N/m');title('Lift Distribution 12,5m/s');
%     
%     
%     
% figure    
%     
%     subplot(3,2,1);
%     plot(y,dD0,'r-o',y,dDpol0,'b');
%     xlabel('y');ylabel('N/m');title('Drag Distribution 0m/s');
%     
%     
%     subplot(3,2,3);
%     plot(y,dD2,'r-o',y,dDpol25,'b');
%     xlabel('y');ylabel('N/m');title('Drag Distribution 2,5m/s');
%     
%     
%     subplot(3,2,5);
%     plot(y,dD50,'r-o',y,dDpol50,'b');
%     xlabel('y');ylabel('N/m');title('Drag Distribution 5m/s');
%     
%     
%     subplot(3,2,2);
%     plot(y,dD75,'r-o',y,dDpol75,'b');
%     xlabel('y');ylabel('N/m');title('Drag Distribution 7,5m/s');
%     
%     
%     subplot(3,2,4);
%     plot(y,dD100,'r-o',y,dDpol100,'b');
%     xlabel('y');ylabel('N/m');title('Drag Distribution 10m/s');
%     
%     
%     subplot(3,2,6);
%     plot(y,dD125,'r-o',y,dDpol125,'b');
%     xlabel('y');ylabel('N/m');title('Drag Distribution 12,5m/s');
  

    

    

    
  



