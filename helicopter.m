clear all;
clc;

%% Data input

% Coses físiques
m = 4500; %[kg]
DL = 350; %[N/m^2]
Mtip = 0.5;
h = 1500; %[m]
vc = 0; % velocitat endavant [m/s]
rroot = 0.1;

% Coses numèriques
nelem = 10;

% Constants
g = 9.81; %[m/s^2]

% Coses òptimes
Cl = 0.615;
Cd = 0.016;
alpha = degtorad(5); %[rad]

nnodes = nelem+1;

%% Pre-cÃ lculs

W = m*g;
A = W/DL; %[m^2]
R = sqrt(A/pi);

% Propietats a l'altitud h (ISA)
Temperature = 288.15-6.5*h/1000; %[K]
rho = 1.225*(Temperature/288.15)^(9810/(6.5*287)-1);

a = sqrt(1.4*287*Temperature);
vtip = a*Mtip; %[m/s]
Omegadisseny = vtip/R; %[rad/s]

%% BEM ideal (MTH)

% Aplicant MTH
vi = sqrt(W/(2*rho*A));
lambdai = vi/(vtip);
lambdac = vc/vtip;

% Discretització
r = linspace(rroot,1,nnodes);
dr = r(2)-r(1);

% Angles
phiideal = atan(lambdai./r);
thetaideal = alpha+phiideal;

% Càlcul de la solidesa
sigmaideal = zeros(1,nnodes);
cideal = zeros(1,nnodes);
nb = 6; % calculat més endavant amb la BEM
for i = 1:nnodes
    sigmaideal(i) = 8*r(i)*(lambdai+lambdac)*lambdai/((r(i)^2+(lambdai+lambdac)^2)*(Cl*cos(phiideal(i))-Cd*sin(phiideal(i))));
    cideal(i) = sigmaideal(i)*pi*R/nb;
    
    % Això no és absolutament necessari
    if sigmaideal(i)<0
        sigmaideal(i) = 0;
    end
    if cideal(i)<0
        cideal(i) = 0;
    end
end

T1=2*rho*pi*R^4*Omegadisseny^2*lambdai^2;
Pi1=T1*vi;
Po1v=0.5*rho*(Omegadisseny*R)^3*pi*R^2*sigmaideal.*r.^3;
Po1=trapz(Po1v);

% figure;
% plot(r,sigmaideal);
% axis([0 1 0 1]);
% xlabel('r')
% ylabel('\sigma')
% title('BEM ideal')
% 
% figure;
% plot(r,cideal);
% xlabel('r')
% ylabel('c (m)')
% title('BEM ideal')
% 
% figure;
% plot(r,thetaideal);
% xlabel('r')
% ylabel('\theta')
% title('BEM ideal')

%% BEM sense pèrdues

for i = 1:nnodes
    if r(i)<=0.7 && r(i+1)>0.7
        l=i;
        break;
    end
end

sigma1 = (sigmaideal(l+1)-sigmaideal(l-1))/(r(l+1)-r(l-1));
sigma0 = sigmaideal(l)-sigma1*r(l);
sigma = zeros(1,nnodes);
for j = 1:nnodes
    sigma(j) = sigma0+sigma1*r(j);
end

theta1 = (thetaideal(l+1)-thetaideal(l-1))/(r(l+1)-r(l-1));

% Per comprovar que l'aproximació lineal de sigma quadra
figure;
plot(r,sigmaideal,r,sigma)
legend('\sigma ideal','\sigma aproximada (r=0.7)')
grid on

lambda = zeros(1,nnodes);
phi = zeros(1,nnodes);
Cl = zeros(1,nnodes);
Cd = zeros(1,nnodes);
thet = zeros(1,nnodes);

Thrust = 0;
thetamin = 0;
thetamax = pi/5;

while abs(Thrust-W)>=1e-2
    
    theta0 = (thetamin+thetamax)/2;
    
    for i = 1:nnodes
        theta = theta0+theta1*r(i);
        thet(i) = theta;
        lambda1 = 0;
        lambda2 = 1;
        phi1 = atan((lambdac+lambda1)/r(i));
        phi2 = atan((lambdac+lambda2)/r(i));
        [Cl1,Cd1] = computeClCd(radtodeg(theta-phi1),0);
        [Cl2,Cd2] = computeClCd(radtodeg(theta-phi2),0);
        F1 = 8*r(i)*lambda1^2-sigma(i)*(r(i)^2+lambda1^2)*(Cl1*cos(phi1)-Cd1*sin(phi1));
        F2 = 8*r(i)*lambda2^2-sigma(i)*(r(i)^2+lambda2^2)*(Cl2*cos(phi2)-Cd2*sin(phi2));
        if F1*F2<0
            F3 = 1;
            while abs(F1-F2)>1e-5
                lambda3 = (lambda1+lambda2)/2;
                phi3 = atan((lambdac+lambda3)/r(i));
                [Cl3,Cd3] = computeClCd(radtodeg(theta-phi3),0);
                F3 = 8*r(i)*lambda3^2-sigma(i)*(r(i)^2+lambda3^2)*(Cl3*cos(phi3)-Cd3*sin(phi3));
                if  F3*F2<0
                    lambda1 = lambda3;
                    F1 = F3;
                elseif  F3*F1<0
                    lambda2 = lambda3;
                    F2 = F3;
                end
            end
        else
            lambda1 = lambda1+0.4;
            lambda2 = lambda2-0.4;
        end
        lambda(i) = lambda3;
        phi(i) = phi3;
        Cl(i) = Cl3;
        Cd(i) = Cd3;
    end
    fprintf('lambdai calculada\n')
    
    Thrust = 0;
    for i = 1:nnodes
        dT = 0.5*rho*vtip^2*(r(i)^2+(lambdac+lambda(i))^2)*(Cl(i)*cos(phi(i))-Cd(i)*sin(phi(i)))*sigma(i)*pi*R^2*dr;
        Thrust = Thrust+dT;
    end
    
    if Thrust-W<0
        thetamin = theta0;
    else
        thetamax = theta0;
    end
    
end

figure;
plot(r,lambda);
xlabel('r');
ylabel('\lambda_{i}');
title('BEM sense pèrdues');

figure;
plot(r,thet);
xlabel('r');
ylabel('\theta (rad)');
title('BEM sense pèrdues');

nb = ceil(sigma(1)*pi*R/0.5);
if nb>=7
    nb = ceil(sigma(1)*pi*R/0.75);
end

c=sigma*pi*R/nb;
figure;
plot(r,c);
xlabel('r');
ylabel('c (m)');
title('BEM sense pèrdues');

Po2v=nb*rho*0.5*(Omegadisseny*R)^2*0.0051*R^2*Omegadisseny*c.*(r.^2).*sqrt(r.^2+lambda.^2);
Po2=trapz(Po2v);

%% BEM+pèrdues+compressibilitat

lambdaPrandtl = zeros(1,nnodes);
FPrandtl = zeros(1,nnodes);
phiPrandtl = zeros(1,nnodes);

while abs(Thrust-W)>=1e-2
    
    theta0 = (thetamin+thetamax)/2;
    
    for i = 1:nnodes
        
        theta = theta0;
        
        % Extrems per fer Bolzano
        lambda1P = 0.1;
        lambda2P = 0.001;
        
        % Càlcul de phi per lambda1 (Prandtl)
        phi21 = atan((lambdac+lambda1P)/r(i));
        phi11 = phi21+1;
        
        while abs(phi11-phi21)>1e-2
            phi11 = phi21;
            F1 = Prandtl(phi11,nb,r(i),rroot);
            phi21 = 0.5*(phi11+atan((lambdac+lambda1P/F1)/r(i)));
        end
        
        % Càlcul de phi per lambda2 (Prandtl)
        phi22 = atan((lambdac+lambda2P)/r(i));
        phi12 = phi22+1;
        
        while abs(phi12-phi22)>1e-2
            phi12 = phi22;
            F2 = Prandtl(phi12,nb,r(i),rroot);
            phi22 = 0.5*(phi12+atan((lambdac+lambda2P/F2)/r(i)));
        end
        
        % Càlcul de les funcions BEM per a lambda1 i lambda2
        [Cl1,Cd1] = computeClCd(radtodeg(theta-phi21),Omegadisseny*r(i)/a);
        [Cl2,Cd2] = computeClCd(radtodeg(theta-phi22),Omegadisseny*r(i)/a);
        Funcio1 = 8*r(i)*(lambda1P+lambdac)*lambda1P-sigma(i)*(r(i)^2+(lambda1P+lambdac)^2)*(Cl1*cos(phi21)-Cd1*sin(phi21));
        Funcio2 = 8*r(i)*(lambda2P+lambdac)*lambda2P-sigma(i)*(r(i)^2+(lambda2P+lambdac)^2)*(Cl2*cos(phi22)-Cd2*sin(phi22));
        
        % Bolzano
        if Funcio1*Funcio2<0
            Funcio3 = 1;
            while abs(Funcio1-Funcio2)>1e-5
                lambda3P = (lambda1P+lambda2P)/2;
                
                % Càlcul de phi per lambda3 (Prandtl)
                phi23 = atan((lambdac+lambda3P)/r(i));
                phi13 = phi23+1;
                while abs(phi13-phi23)>1e-2
                    phi13 = phi23;
                    F3 = Prandtl(phi13,nb,r(i),rroot);
                    phi23 = 0.5*(phi13+atan((lambdac+lambda3P/F3)/r(i)));
                end
                
                % Càlcul de la funció BEM per lambda3
                [Cl3,Cd3] = computeClCd(radtodeg(theta-phi23),Omegadisseny*r(i)/a);
                Funcio3 = 8*r(i)*(lambda3P+lambdac)*lambda3P-sigma(i)*(r(i)^2+(lambda3P+lambdac)^2)*(Cl3*cos(phi23)-Cd3*sin(phi23));
                
                % Comparació per fer Bolzano
                if  Funcio3*Funcio2<0
                    lambda1P = lambda3P;
                    Funcio1 = Funcio3;
                elseif  Funcio3*Funcio1<0
                    lambda2P = lambda3P;
                    Funcio2 = Funcio3;
                end
            end
        else
            % no assigna valors als punts en què Bolzano no convergeix
            lambda3P = NaN;
            F3 = NaN;
            phi23 = NaN;
        end
        
        lambdaPrandtl(i) = lambda3P;
        FPrandtl(i) = F3;
        phiPrandtl(i) = phi23;
        ClPrandtl(i) = Cl3;
        CdPrandtl(i) = Cd3;
    end
    
    fprintf('lambdai calculada\n')
    
    Thrust = 0;
    for i = 1:nnodes
        dT = 0.5*rho*vtip^2*(r(i)^2+(lambdac+lambda(i))^2)*(Cl(i)*cos(phi(i))-Cd(i)*sin(phi(i)))*sigma(i)*pi*R^2*dr;
        Thrust = Thrust+dT;
    end
    
    if Thrust-W<0
        thetamin = theta0;
    else
        thetamax = theta0;
    end
    
end

% figure;
% plot(r,lambda,r,lambdaPrandtl);
% title('BEM+pèrdues');
% legend('\lambda sense pèrdues','\lambda amb pèrdues');
% 
% figure;
% plot(r,FPrandtl);
% xlabel('r');
% ylabel('F');