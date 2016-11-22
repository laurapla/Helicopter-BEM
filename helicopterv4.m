clear all;
clc;

%% Data input

% Coses físiques
m = 4500; %[kg]
DL = 350; %[N/m^2]
Mtip = 0.5;
h = 1500; %[m]
rroot = 0.1;
n=0;
nn=0;

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
Vc = [0 2.5 5 7.5 10 12.5];

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
lambdai1 = vi/(vtip);
lambdac = Vc / vtip;

% Discretització
r = linspace(rroot,1,nnodes);
dr = r(end)/nelem;

% Càlcul de angles i solidesa
% Cada columna de la matriu conte la variable calculada per a una velocitat de climbing
sigmaideal = zeros(nnodes,length(Vc));
cideal = zeros(nnodes,length(Vc));
nb = 6; % calculat més endavant amb la BEM

for j=1:length(Vc)  % Calcula la solides i angles per tots els valors de Vc
    for i = 1:nnodes
        
        phiideal(i,j) = atan((lambdai1+lambdac(j))/r(i));
        thetaideal(i,j) = alpha0+phiideal(i,j);
        
        sigmaideal(i,j) = 8*r(i)*lambdai1*(lambdai1+lambdac(j))/((r(i)^2+(lambdai1+lambdac(j))^2)*(Cl*cos(phiideal(i,j))-Cd*sin(phiideal(i,j))));
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
       Pi(i,j) = W*vtip*(lambdai1 + lambdac(j));
   end
end

for i=1:length(Vc)
   Pi_total(i) = trapz(r,Pi(:,i));
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
   
    for jm = 1:nnodes
      sigma(jm,j) = sigma0(j)+sigma1(j)*r(jm);
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

ewlambdai = 0;  
thetamax = 20;
thetamin = -5;
Thrust = zeros(1,length(Vc));         % inicialitzacio del Thrust        
Pi = zeros(1,length(Vc));    % inicialitzaci� de la pot�ncia indu�da [W]
P0 = zeros(1,length(Vc));  % inicialitzci� de la pot�ncia par�sita [W]c


for j=1:length(Vc)
    
    while abs(Thrust(j)-W)>er
        Thrust(j)-W
       
        theta0it = 0.5*(thetamax+thetamin);   % calcul del punt mig de l'interval
        
        % Camp de velocitats induides
        
        
        for i=1:nnodes
            
            i
            lambdait1 = lambdai1;
            lambdait2 = lambdai1;
            
            theta(i,j) = theta0it*(pi/180)+theta1(j)*r(i);  thetadeg(i,j) = theta(i,j)*(180/pi);  % angle theta [deg]
            F = 9;
            
            while abs(F)> er2
                
                phiit= atan((lambdait1+lambdac(j))/r(i)); phiitdeg=phiit*(180/pi);
                alphait= thetadeg(i,j)-phiitdeg;
                Mit = vtip * sqrt( ( (lambdait1+lambdac(j))^2 + r(i)^2 ) )/a;
                [clit,cdit] = computeClCd(alphait, Mit);
                Kit = clit*cos(phiit)-cdit*sin(phiit);
                AA=8*r(i); BB=lambdac(j); CC=r(i)^2; EE = sigma(i,j)* Kit;
                
                F1 = AA*(BB*+lambdait1)*lambdait1-(CC*+(lambdait1+BB)^2)*EE;
                
                
                phiit= atan((lambdait2+lambdac(j))/r(i)); phiitdeg=phiit*(180/pi);
                alphait= thetadeg(i,j)-phiitdeg;
                Mit = vtip * sqrt( ( (lambdait2+lambdac(j))^2 + r(i)^2 ) )/a;
                [clit,cdit] = computeClCd(alphait, Mit);
                Kit = clit*cos(phiit)-cdit*sin(phiit);
                AA=8*r(i); BB=lambdac(j); CC=r(i)^2; EE = sigma(i,j)* Kit;
                
                F2 = AA*(BB*+lambdait2)*lambdait2-(CC*+(lambdait2+BB)^2)*EE;
                
                
                if F1*F2 <=0
                    
                    lambdait = 0.5*(lambdait1+lambdait2);
                    
                    phiit= atan((lambdait+lambdac(j))/r(i)); phiitdeg=phiit*(180/pi);
                    alphait= thetadeg(i,j)-phiitdeg;
                    Mit = vtip * sqrt( ( (lambdait+lambdac(j))^2 + r(i)^2 ) )/a;
                    [clit,cdit] = computeClCd(alphait, Mit);
                    Kit = clit*cos(phiit)-cdit*sin(phiit);
                    AA=8*r(i); BB=lambdac(j); CC=r(i)^2; EE = sigma(i,j)* Kit;
                    
                    F = AA*(BB*+lambdait)*lambdait-(CC*+(lambdait+BB)^2)*EE;
                    
                    if abs(F)<= er2
                        lambda(i,j) = lambdait;
                    else
                        if F*F1<0
                            lambdait2 = lambdait; 
                        else
                            lambdait1 = lambdait; 
                        end
                    end
                    
                else
                    
                    lambdait2 = lambdait2+0.0002;
                    lambdait1 = lambdait1-0.0002;
                end
                
            end
            
            phi(i,j) = atan((lambda(i,j)+lambdac(j))/r(i));  phideg(i,j) = phi(i,j)*(180/pi);            % angle de la velocitat incident [deg]
            alpha(i,j) = theta(i,j)-phi(i,j);  alphadeg(i,j) = alpha(i,j)*(180/pi);                  % angle d'atac [deg]
            M(i,j) = vtip*sqrt( (lambda(i,j)+lambdac(j))^2 + r(i)^2 )/a;                                    % numero de Mach a l'element de pala
            
            [cl(i,j),cd(i,j)] = computeClCd(alphadeg(i,j), M(i,j));               % funci� que d�na Cl i Cd en funci� d'alpha i Mach
            
            K2 = cl(i,j)*cos(phi(i,j))-cd(i,j)*sin(phi(i,j));  % coeficient utilitzat per partir l'expressi� del dT
            
            dT(j) = nb*0.5*rho*c(i,j)*vtip^2*R*(r(i)^2+(lambda(i,j)+lambdac(j))^2)*K2*dr;   % Calucl diferencial de Thrust
            Thrust(j) = Thrust(j)+dT(j);            % Thrust de l'helicopter [N]
            
            dPi(j) = vtip*(lambda(i,j) + lambdac(j))*dT(j);   % Calcul diferencial de potencia induida
            Pi(j) = Pi(j)+dPi(j);                          % Potencia induida
        end
        
        if abs(Thrust(j)-W)> er
           
            if Thrust(j)-W>=0
               thetamax = theta0it;
            else
            thetamin = theta0it;
            
            end
        else
            theta0(j)=theta0it;
        end
    end
end



% Po2v=nb*rho*0.5*(Omegadisseny*R)^2*0.0051*R^2*Omegadisseny*c.*(r.^2).*sqrt(r.^2+lambda.^2);
% Po2=trapz(Po2v);
% 
% 
% figure
% 
% grid on
% plot(Vc, Pi_total);
% 
% 
% figure;
% plot(r,lambda);
% xlabel('r');
% ylabel('\lambda_{i}');
% title('BEM sense pèrdues');
