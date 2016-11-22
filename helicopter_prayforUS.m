clear all;
clc;

%% Data input

% Coses f√≠siques
m = 4500; %[kg]
DL = 350; %[N/m^2]
Mtip = 0.5;
h = 1500; %[m]

% Coses num√®riques
nelem = 20;

% Constants
g = 9.81; %[m/s^2]

% Coses √≤ptimes 
Cl = 1.0248;
Cd = 0.01077;
alpha0 = degtorad(9); %[rad]

nnodes = nelem+1;

% Velocitats d'ascens
Vc = [0 2.5 5 7.5 10 12.5]; nVc=length(Vc);

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
lambdac = Vc / vtip;

% Discretitzaci√≥
r = linspace(0.1,1,nnodes);
dr = r(end)/nelem;

% C√†lcul de angles i solidesa
% Cada columna de la matriu conte la variable calculada per a una velocitat de climbing
sigmaideal = zeros(nnodes,nVc);
cideal = zeros(nnodes,nVc);
nb = 6; % calculat m√©s endavant amb la BEM

for j=1:nVc  % Calcula la solides i angles per tots els valors de Vc
    for i = 1:nnodes
        
        phiideal(i,j) = atan((lambdai+lambdac(j))/r(i));
        thetaideal(i,j) = alpha0+phiideal(i,j);
        
        sigmaideal(i,j) = 8*r(i)*lambdai*(lambdai+lambdac(j))/((r(i)^2+(lambdai+lambdac(j))^2)*(Cl*cos(phiideal(i,j))-Cd*sin(phiideal(i,j))));
        cideal(i,j) = sigmaideal(i,j)*pi*R/nb;
    
        % Aix√≤ no √©s absolutament necessari
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
   for j=1:nVc 
       Pi(i,j) = W*vtip*(lambdai + lambdac(j));
   end
end

for i=1:nVc
   Pi_total(i) = trapz(r.*R,Pi(:,i));
end

%Calcul de la Potencia Parasita

figure;
plot(r,sigmaideal);
axis([0 1 0 1]);
xlabel('r')
ylabel('\sigma')
title('SIGMA - BEM ideal')

figure;
plot(r,cideal);
xlabel('r')
ylabel('c (m)')
title('CORDA - BEM ideal')

figure;
plot(r,thetaideal);
xlabel('r')
ylabel('\theta')
title('THETA - BEM ideal')

fprintf('MTH computed.\n')

%% BEM sense p√®rdues
for i = 1:nnodes
    if r(i)<=0.7 && r(i+1)>0.7
        l=i;    %% Node al qual r=0.7
        break;
    end
end

%Calcul solidesa i torsio per cada Vc
sigma = zeros(nnodes,nVc);
for j=1:length(Vc)
    sigma1(j) = (sigmaideal(l+1,j)-sigmaideal(l-1,j))/(r(l+1)-r(l-1));
    sigma0(j) = sigmaideal(l,j)-sigma1(j)*r(l);
    
    theta1(j) = (thetaideal(l+1,j)-thetaideal(l-1,j))/(r(l+1)-r(l-1));
%   theta0(j) = alpha0+phiideal(l,j)*-theta1(j)*r(l);
    theta0(j) = 0; % Unknown of the problem
    
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

%Calcul distribucio lift


figure;
plot(r,c);
xlabel('r')
ylabel('c (m)')
title('CORDA - BEM no ideal')

%% Computation of thetaC for each of the Vc
prec_thrust = 1e-4;  
prec_lambda=1e-4;       % error maxim permes
error_lambda = 100; 

thetaT = zeros(nnodes,nVc);
lambda = zeros(nnodes,nVc); phideg=lambda; alpha=lambda;

fprintf('Starting BEM solution\n')

for j=1:length(Vc)

    Thrust = 0;         % inicialitzacio del Thrust
    lambdNEW = 0;  
    thetaC_min = 2;     % angle de pas colectiu minim [deg] (aquell que fa que W>T)
    thetaC_max = 5;     % angle de pas colectiu maxim [deg] (aquell que fa que W<T)

    while abs(Thrust-W)>prec_thrust

        thetaC = 0.5*(thetaC_max+thetaC_min);   % calcul del punt mig de l'interval 
        theta_r=theta(:,j)*180/pi+thetaC;       % Vector de distribuciÛ de colectiu
        
        % Camp de velocitats induides
        for i=1:nnodes
            
%             % QuÈ pinta tiene F, y lambda?
%             lambda_test=-20:0.1:20; F=zeros(1,length(lambda_test));
%             for t=1:length(lambda_test)
%                 F(t)=evaluateBEM(lambda_test(t),lambdac(j),r(i),theta_r(i),sigma(i,j));
%             end
%             plot(lambda_test,F)
            
            lambda1 = -1;          % inicialitzem de lambda a la del rotor ideal
            lambda2 = 1;
            
            F1=evaluateBEM(lambda1,lambdac(j),r(i),theta_r(i),sigma(i,j));
            F2=evaluateBEM(lambda2,lambdac(j),r(i),theta_r(i),sigma(i,j));
            
            if F1*F2<0
                F3 = 1;
                
                for n=1:1000; % Iterem per cada posicio i
                lambda3 = (lambda1+lambda2)/2;
                F3 = evaluateBEM(lambda3,lambdac(j),r(i),theta_r(i),sigma(i,j));
                    if F3*F2<0
                        lambda1=lambda3;
                        F1=evaluateBEM(lambda1,lambdac(j),r(i),theta_r(i),sigma(i,j));
                    elseif F3*F1<0
                        lambda2=lambda3;
                        F2=evaluateBEM(lambda2,lambdac(j),r(i),theta_r(i),sigma(i,j));
                    end
                
                % EvaluaciÛ de la nova lamda calculada
                if abs(F1-F2)<prec_lambda
                    break
                elseif n==10000
                    fprintf('\n Position r(%g)=%g unconverged\n',i,r(i))
                end
%                 if imag(lambdNEW)~=0 % Recagada
%                     fprintf('Lambda de la it=%g imaginaria tete. Go home you are drunk\n',n)
%                     break
%                 elseif abs(error_lambda)>prec_lambda % Tornem a iterar
%                     % lambdait = 0.5*(lambdNEW+lambdait);   
%                     lambdait = lambdNEW;
%                 elseif n==1000
%                     fprintf('\n Position i=%g unconverged\n',i)
%                 else % converged! Enregistrem lambdaNEW
%                     lambda(i,j)=lambdNEW;
%                     phideg(i,j)=phi_n;
%                     alpha(i,j)=alpha_n;
%                 end
                end % End for de n iteracions per lambda(r(i))
            end
            lambda(i,j)=lambda3;
        end % End de tots els nodes
         
        % VerificaciÛ de correcte c‡lcul de velocitats induides
        plot(r,lambda(:,j))
        
        
    
    
    % Calcul del Thrust i les potÔøΩncies induÔøΩdes i parÔøΩsites
    Thrust = 0; % inicialitzaciÔøΩ de l'empenta [N]
    Pi = 0;     % inicialitzaciÔøΩ de la potÔøΩncia induÔøΩda [W]
    P0 = 0;     % inicialitzaciÔøΩ de la potÔøΩncia parÔøΩsita [W]
    
    % Guardem la soluciÛ
    thetaT(i,j) = theta(i,j)*(180/pi)+thetaC;    % angle theta [deg] 
    end
end

%% Tornem a repetir-ho tot, pk podem
 for j=1:length(Vc)
    for i=1:nnodes

        phi_n(i,j) = atan(lambda(i,j)/r(i)); phideg(i,j) = phi_n(i,j)*(180/pi);            % angle de la velocitat incident [deg]
        alphadeg (i,j) = thetaT(i,j)-phideg(i,j);  alpha(i,j) = alphadeg(i,j)*(pi/180);                  % angle d'atac [deg]
        M = vtip*sqrt((lambda(i,j)+lambdac(j))^2+r(i)^2)/a; % numero de Mach a l'element de pala

        [cl(i,j),cd(i,j)] = computeClCd(alphadeg(i,j), M(i,j));               % funciÔøΩ que dÔøΩna Cl i Cd en funciÔøΩ d'alpha i Mach

        K2 = cl(i,j)*cos(phi_n(i,j))-cd(i,j)*sin(phi_n(i,j));  % coeficient utilitzat per partir l'expressiÔøΩ del dT

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


 

Po2v=nb*rho*0.5*(Omegadisseny*R)^2*0.0051*R^2*Omegadisseny*c.*(r.^2).*sqrt(r.^2+lambda.^2);
Po2=trapz(Po2v);


% figure
% 
% grid on
% plot(Vc, Pi_total);


% figure;
% plot(r,lambda);
% xlabel('r');
% ylabel('\lambda_{i}');
% title('BEM sense p√®rdues');

