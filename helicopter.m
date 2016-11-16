clear all;
clc;

%% Data input

% Coses físiques
m = 4500; %[kg]
DL = 350; %[N/m^2]
Mtip = 0.5;
h = 1500; %[m]

% Coses numèriques
nelem = 20;

% Constants
g = 9.81; %[m/s^2]

% Coses òptimes 
Cl = 0.615;
Cd = 0.016;
alpha = degtorad(5); %[rad]

nnodes = nelem+1;

%% Pre-càlculs

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

% Discretització
r = linspace(0,1,nnodes);

% Angles
phiideal = atan(lambdai./r);
thetaideal = alpha+phiideal;

% Càlcul de la solidesa
sigmaideal = zeros(1,nnodes);
cideal = zeros(1,nnodes);
nb = 7;
for i = 1:nnodes
    sigmaideal(i) = 8*r(i)*lambdai^2/((r(i)^2+lambdai^2)*(Cl*cos(phiideal(i))-Cd*sin(phiideal(i))));
    cideal(i) = sigmaideal(i)*pi*R/nb;
    
    % Això no és absolutament necessari
    if sigmaideal(i)<0
        sigmaideal(i) = 0;
    end
    if cideal(i)<0
        cideal(i) = 0;
    end
end

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
        break;
    end
end

sigma1 = (sigmaideal(i+1)-sigmaideal(i-1))/(r(i+1)-r(i-1));
sigma0 = sigmaideal(i)-sigma1*r(i);
sigma = zeros(1,nnodes);
for jm = 1:nnodes
    sigma(jm) = sigma0+sigma1*r(jm);
end

theta1 = (thetaideal(i+1)-thetaideal(i-1))/(r(i+1)-r(i-1));

% Per comprovar que l'aproximació lineal de sigma quadra
figure;
plot(r,sigmaideal,r,sigma)
legend('\sigma ideal','\sigma aproximada (r=0.7)')
grid on

lambda = zeros(1,nnodes);

for i = 1:nnodes
    lambda1 = 1;
    lambda2 = -0.001;
    F1 = 8*r(i)*lambda1^2-sigma(i)*(r(i)^2+lambda1^2)*(Cl*cos(atan(lambda1/r(i)))-Cd*sin(atan(lambda1/r(i))));
    F2 = 8*r(i)*lambda2^2-sigma(i)*(r(i)^2+lambda2^2)*(Cl*cos(atan(lambda2/r(i)))-Cd*sin(atan(lambda2/r(i))));
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
title('BEM sense pèrdues');