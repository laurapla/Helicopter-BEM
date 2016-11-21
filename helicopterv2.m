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
alpha0 = degtorad(5); %[rad]

nnodes = nelem+1;

% Velocitats d'ascens
Vc = [2.5 5 7.5 10 12.5];


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
lambdac = Vc / vtip;

% Discretització
r = linspace(0,1,nnodes);
dr = r(end)/nelem;

% Càlcul de angles i solidesa
% Cada columna de la matriu conte la variable calculada per a una velocitat de climbing
sigmaideal = zeros(nnodes,length(Vc));
cideal = zeros(nnodes,length(Vc));
nb = 6; % calculat més endavant amb la BEM

for j=1:length(Vc)  % Calcula la solides i angles per tots els valors de Vc
    for i = 1:nnodes
        
        phiideal(i,j) = atan((lambdai+lambdac(j))/r(i));
        thetaideal(i,j) = alpha0+phiideal(i,j);
        
        sigmaideal(i,j) = 8*r(i)*lambdai*(lambdai+lambdac(j))/((r(i)^2+(lambdai+lambdac(j))^2)*(Cl*cos(phiideal(i,j))-Cd*sin(phiideal(i,j))));
        cideal(i,j) = sigmaideal(i,j)*pi*R/nb;
    
        % Això no és absolutament necessari
        if sigmaideal(i,j)<0
        sigmaideal(i,j) = 0;
        end
        if cideal(i,j)<0
        cideal(i,j) = 0;
        end
    end
end

%Calcul de la Potencia induida

for i=1:nnodes 
   for j=1:length(Vc) 
       Pi(i,j) = W*vtip*(lambdai + lambdac(j));
   end
end

for i=1:length(Vc)
   Pi_total(i) = trapz(Pi(:,i));
end

%Calcul de la Potencia Parasita



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
        l=i;    %% Node al qual r=0.7
        break;
    end
end

%Calcul solidesa i torsio
sigma = zeros(nnodes,length(Vc));
for j=1:length(Vc)
    sigma1(j) = (sigmaideal(l+1,j)-sigmaideal(l-1,j))/(r(l+1)-r(l-1));
    sigma0(j) = sigmaideal(l,j)-sigma1(j)*r(l);
    
    theta1(j) = (thetaideal(l+1,j)-thetaideal(l-1,j))/(r(l+1)-r(l-1));
    theta0 (j) = alpha0+phiideal(l,j)*-theta1(j)*r(l);
    
    for jm = 1:nnodes
      sigma(jm,j) = sigma0(j)+sigma1(j)*r(jm);
      theta(jm,j) = theta0(j)+theta1(j)*r(jm);
    end
end

nb = ceil(sigma(1,1)*pi*R/0.5);
if nb>=7
    nb = ceil(sigma(1,1)*pi*R/0.75);
end

c=sigma*pi*R/nb;

% Calcul de lambda

er = 1e-4;  er2=1e-8;       % error maxim permes
error = 100; 

lambdNEW = 0;  
thetaC_min = 2;     % angle de pas colectiu minim [deg] (aquell que fa que W>T)
thetaC_max = 5;     % angle de pas colectiu maxim [deg] (aquell que fa que W<T)

Thrust = 0;         % inicialitzacio del Thrust         

while abs(Thrust-W)>er
    
    thetaC = 0.5*(thetaC_max+thetaC_min);   % calcul del punt mig de l'interval 
    
        % Camp de velocitats induides
 for j=1:length(Vc)
    for i=1:nnodes

        lambdait = lambdai;          % inicialitzaci� de lambda a la de Hovering [-]
        thetaT(i,j) = theta(i,j)*(180/pi)+thetaC;    % angle theta [deg]  
        error = 2;
            
        while abs(error)> er2
    
            phideg(i,j) = atan((lambdait+lambdac(j))/r(i))*(180/pi);               % angle de la velocitat incident [deg]
            alpha(i,j) = thetaT(i,j)-phideg(i,j);                              % angle d'atac [deg]
            M(i,j) = vtip * sqrt( ( (lambdait+lambdac(j))^2 + r(i)^2 ) )/a;         % numero de Mach a l'element de pala
            [cl(i,j),cd(i,j)] = computeClCd(alpha(i,j), M(i,j));                       % coeficients aerodinamics de l'element de pala
        
            K(i,j) = cl(i,j)*cos(phiideal(i,j))-Cd*sin(phiideal(i,j));
            
            AA=8*r(i); BB=lambdac(j); CC=r(i)^2; EE = sigma(i,j)* K(i,j);
           
            lambdNEW = (-BB*(AA-2*EE)+ sqrt( BB^2 * (AA-2*EE)^2 - 4*(AA-EE)*EE*(CC-BB^2)))/(2*(AA-EE)); 
            
            error = lambdait-lambdNEW;
        
            if abs(error)>er2
                lambdait = lambdait+0.5*(lambdNEW-lambdait);
            else
                lambda(i,j) = 0.5*(lambdait+lambdNEW);
       
            end
        end
   
    end
end

        % Calcul del Thrust i les pot�ncies indu�des i par�sites
    Thrust = 0; % inicialitzaci� de l'empenta [N]
    Pi = 0;     % inicialitzaci� de la pot�ncia indu�da [W]
    P0 = 0;     % inicialitzaci� de la pot�ncia par�sita [W]
    
 for j=1:length(Vc)
    Thrust(j)=0; 
    Pi(j)=0;
    for i=1:nnodes

        phi(i,j) = atan(lambda(i,j)/r(i)); phideg(i,j) = phi(i,j)*(180/pi);            % angle de la velocitat incident [deg]
        alphadeg (i,j) = thetaT(i,j)-phideg(i,j);  alpha(i,j) = alphadeg(i,j)*(pi/180);                  % angle d'atac [deg]
       
        [cl(i,j),cd(i,j)] = computeClCd(alphadeg(i,j), M(i,j));               % funci� que d�na Cl i Cd en funci� d'alpha i Mach

        K2 = cl(i,j)*cos(phi(i,j))-cd(i,j)*sin(phi(i,j));  % coeficient utilitzat per partir l'expressi� del dT

        dT(j) = nb*0.5*rho*c(i,j)*vtip^2*R*(r(i)^2+(lambda(i,j)+lambdac(j))^2)*K2*dr;   % Calucl diferencial de Thrust
        Thrust(j) = Thrust(j)+dT(j);            % Thrust de l'helicopter [N]
        
        dPi(j) = vtip*(lambda(i,j) + lambdac(j))*dT(j);   % Calcul diferencial de potencia induida
        Pi(j) = Pi(j)+dPi(j);                                 % Potencia induida
    end
   
    if Thrust-W>=0
        thetaC_max = thetaC;        
    else
        thetaC_min = thetaC;        
            
    end
 end
end

 

Po2v=nb*rho*0.5*(Omegadisseny*R)^2*0.0051*R^2*Omegadisseny*c.*(r.^2).*sqrt(r.^2+lambda.^2);
Po2=trapz(Po2v);


figure

grid on
plot(Vc, Pi_total);


figure;
plot(r,lambda);
xlabel('r');
ylabel('\lambda_{i}');
title('BEM sense pèrdues');

