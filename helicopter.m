clear all;
clc;

%% Data input

% Coses f√≠siques
m = 4500; %[kg]
DL = 350; %[N/m^2]
Mtip = 0.5;
h = 1500; %[m]
vc = 0; % velocitat endavant [m/s]
rroot = 0.1;

% Coses num√®riques
nelem = 100;

% Constants
g = 9.81; %[m/s^2]

% Coses √≤ptimes 
Cl = 0.615;
Cd = 0.016;
alpha = degtorad(5); %[rad]

nnodes = nelem+1;

%% Pre-c√†lculs

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

% Discretitzaci√≥
r = linspace(rroot,1,nnodes);

% Angles
phiideal = atan(lambdai./r);
thetaideal = alpha+phiideal;

% C√†lcul de la solidesa
sigmaideal = zeros(1,nnodes);
cideal = zeros(1,nnodes);
nb = 6; % calculat m√©s endavant amb la BEM
for i = 1:nnodes
    sigmaideal(i) = 8*r(i)*lambdai^2/((r(i)^2+lambdai^2)*(Cl*cos(phiideal(i))-Cd*sin(phiideal(i))));
    cideal(i) = sigmaideal(i)*pi*R/nb;
    
    % Aix√≤ no √©s absolutament necessari
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

%% BEM sense p√®rdues

for i = 1:nnodes
    if r(i)<=0.7 && r(i+1)>0.7
        l=i;
        break;
    end
end

sigma1 = (sigmaideal(l+1)-sigmaideal(l-1))/(r(l+1)-r(l-1));
sigma0 = sigmaideal(l)-sigma1*r(l);
sigma = zeros(1,nnodes);
for jm = 1:nnodes
    sigma(jm) = sigma0+sigma1*r(jm);
end

theta1 = (thetaideal(l+1)-thetaideal(l-1))/(r(l+1)-r(l-1));

% Per comprovar que l'aproximaci√≥ lineal de sigma quadra
figure;
plot(r,sigmaideal,r,sigma)
legend('\sigma ideal','\sigma aproximada (r=0.7)')
grid on

lambda = zeros(1,nnodes);

for i = 1:nnodes
    lambda1 = 0.1;
    lambda2 = 0.01;
    phi1 = atan((lambdac+lambda1)/r(i));
    phi2 = atan((lambdac+lambda2)/r(i));
    F1 = 8*r(i)*lambda1^2-sigma(i)*(r(i)^2+lambda1^2)*(Cl*cos(phi1)-Cd*sin(phi1));
    F2 = 8*r(i)*lambda2^2-sigma(i)*(r(i)^2+lambda2^2)*(Cl*cos(phi2)-Cd*sin(phi2));
    if F1*F2<0
        F3 = 1;
        while abs(F1-F2)>1e-5
            lambda3 = (lambda1+lambda2)/2;
            F3 = 8*r(i)*lambda3^2-sigma(i)*(r(i)^2+lambda3^2)*(Cl*cos(atan(lambda3/r(i)))-Cd*sin(atan(lambda3/r(i))));
            if  F3*F2<0
                lambda1 = lambda3;
                F1 = 8*r(i)*lambda1^2-sigma(i)*(r(i)^2+lambda1^2)*(Cl*cos(atan(lambda1/r(i)))-Cd*sin(atan(lambda1/r(i))));
            elseif  F3*F1<0
                lambda2 = lambda3;
                F2 = 8*r(i)*lambda2^2-sigma(i)*(r(i)^2+lambda2^2)*(Cl*cos(atan(lambda2/r(i)))-Cd*sin(atan(lambda2/r(i))));
            end
        end
    end
    lambda(i) = lambda3;
end

figure;
plot(r,lambda);
xlabel('r');
ylabel('\lambda_{i}');
title('BEM sense p√®rdues');

nb = ceil(sigma(1)*pi*R/0.5);
if nb>=7
    nb = ceil(sigma(1)*pi*R/0.75);
end

c=sigma*pi*R/nb;
Po2v=nb*rho*0.5*(Omegadisseny*R)^2*0.0051*R^2*Omegadisseny*c.*(r.^2).*sqrt(r.^2+lambda.^2);
Po2=trapz(Po2v);

%% BEM+pËrdues+compressibilitat

lambdaPrandtl = zeros(1,nnodes);
FPrandtl = zeros(1,nnodes);
phiPrandtl = zeros(1,nnodes);

for i = 1:nnodes
    
    % C‡lcul dels coeficients
    [Cl,Cd] = computeClCd(5,Omegadisseny*r(i)/a);
    
    % Extrems per fer Bolzano
    lambda1P = 0.1;
    lambda2P = 0.01;
    
    % C‡lcul de phi per lambda1 (Prandtl)
    phi21 = atan((lambdac+lambda1P)/r(i));
    phi11 = phi21+1;
    
    while abs(phi11-phi21)>1e-2
        phi11 = phi21;
        F1 = Prandtl(phi11,nb,r(i),rroot);
        phi21 = 0.5*(phi11+atan((lambdac+lambda1P/F1)/r(i)));
    end
    
    % C‡lcul de phi per lambda2 (Prandtl)
    phi22 = atan((lambdac+lambda2P)/r(i));
    phi12 = phi22+1;
    
    while abs(phi12-phi22)>1e-2
        phi12 = phi22;
        F2 = Prandtl(phi12,nb,r(i),rroot);
        phi22 = 0.5*(phi12+atan((lambdac+lambda2P/F2)/r(i)));
    end
    
    % C‡lcul de les funcions BEM per a lambda1 i lambda2
    Funcio1 = 8*r(i)*lambda1P^2-sigma(i)*(r(i)^2+lambda1P^2)*(Cl*cos(phi21)-Cd*sin(phi21));
    Funcio2 = 8*r(i)*lambda2P^2-sigma(i)*(r(i)^2+lambda2P^2)*(Cl*cos(phi22)-Cd*sin(phi22));
    
    % Bolzano
    if Funcio1*Funcio2<0
        Funcio3 = 1;
        while abs(Funcio1-Funcio2)>1e-5
            lambda3P = (lambda1P+lambda2P)/2;
            
            % C‡lcul de phi per lambda3 (Prandtl)
            phi23 = atan((lambdac+lambda3P)/r(i));
            phi13 = phi23+1;
            while abs(phi13-phi23)>1e-2
                phi13 = phi23;
                F3 = Prandtl(phi13,nb,r(i),rroot);
                phi23 = 0.5*(phi13+atan((lambdac+lambda3P/F3)/r(i)));
            end
            
            % C‡lcul de la funciÛ BEM per lambda3
            Funcio3 = 8*r(i)*lambda3P^2-sigma(i)*(r(i)^2+lambda3P^2)*(Cl*cos(phi23)-Cd*sin(phi23));
            
            % ComparaciÛ per fer Bolzano
            if  Funcio3*Funcio2<0
                lambda1P = lambda3P;
                Funcio1 = Funcio3;
            elseif  Funcio3*Funcio1<0
                lambda2P = lambda3P;
                Funcio2 = Funcio3;
            end
        end
    else
        % no assigna valors als punts en quË Bolzano no convergeix
        lambda3P = NaN;
        F3 = NaN;
        phi23 = NaN;
    end
    
    lambdaPrandtl(i) = lambda3P;
    FPrandtl(i) = F3;
    phiPrandtl(i) = phi23;
end

figure;
plot(r,lambda,r,lambdaPrandtl);
title('BEM+pËrdues');
legend('\lambda sense pËrdues','\lambda amb pËrdues');

figure;
plot(r,FPrandtl);
xlabel('r');
ylabel('F');