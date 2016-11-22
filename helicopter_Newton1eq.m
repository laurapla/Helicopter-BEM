clear all;
clc;

%% Data input

% Coses fÃ­siques
m = 4500; %[kg]
DL = 350; %[N/m^2]
Mtip = 0.5;
h = 1500; %[m]
rroot = 0.1;

% Coses numÃ¨riques
nelem = 20;

% Constants
g = 9.81; %[m/s^2]

% Coses Ã²ptimes 
Cl = 1.0248;
Cd = 0.01077;
alpha0 = degtorad(9); %[rad]

nnodes = nelem+1;

% Velocitats d'ascens
Vc = [0 2.5 5 7.5 10 12.5]; nVc=length(Vc);

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
lambdac = Vc / vtip;

% DiscretitzaciÃ³
r = linspace(rroot,1,nnodes);
dr = r(end)-r(end-1);

% CÃ lcul de angles i solidesa
% Cada columna de la matriu conte la variable calculada per a una velocitat de climbing
sigmaideal = zeros(nnodes,nVc);
cideal = zeros(nnodes,nVc);
nb = 6; % calculat mÃ©s endavant amb la BEM

for j=1:nVc  % Calcula la solides i angles per tots els valors de Vc
    for i = 1:nnodes
        
        phiideal(i,j) = atan((lambdai+lambdac(j))/r(i));
        thetaideal(i,j) = alpha0+phiideal(i,j);
        
        sigmaideal(i,j) = 8*r(i)*lambdai*(lambdai+lambdac(j))/((r(i)^2+(lambdai+lambdac(j))^2)*(Cl*cos(phiideal(i,j))-Cd*sin(phiideal(i,j))));
        cideal(i,j) = sigmaideal(i,j)*pi*R/nb;
    
        % AixÃ² no Ã©s absolutament necessari
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

% figure;
% plot(r,sigmaideal);
% axis([0 1 0 1]);
% xlabel('r')
% ylabel('\sigma')
% title('SIGMA - BEM ideal')
% 
% figure;
% plot(r,cideal);
% xlabel('r')
% ylabel('c (m)')
% title('CORDA - BEM ideal')
% 
% figure;
% plot(r,thetaideal);
% xlabel('r')
% ylabel('\theta')
% title('THETA - BEM ideal')

fprintf('MTH computed.\n')

%% BEM sense pÃ¨rdues
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

%% NEWTON METHOD FOR Computation of thetaC for each of the Vc
% Numerical constants
prec_thrust = 1e-4;  
prec_lambda=1e-8;       % error maxim permes
error_lambda = 100; 
nmaxIT=100;
nmaxITtheta=15;

% Preallocation of variables
thetaT = zeros(nnodes,nVc);
lambda = zeros(nnodes,nVc); phideg=lambda; alpha=lambda;
CL = zeros(nnodes,1); 
CD = zeros(nnodes,1); 
lc=lambdac(j);  rf=r(i); sigf=sigma(i,j); cl=0.5; cd=cl; dcl=cl; dcd=cl;

% Functions arrays
f=@(li,fi) 8*rf*(lc+li)*li-...
    (rf^2+(lc+li)^2)*sigf*(cl*cosd(fi)-cd*sind(fi));

% Jacobian of the functions
Df=@(li,fi) 8*lc*rf+16*li*rf-2*(lc+li)*sigf*(cl*cosd(fi)-cd*sind(fi))+...
    sigf*rf^2*(dcl*cosd(fi)+cl*sind(fi)-dcd*sind(fi)+cd*cosd(fi)) ;

fprintf('Starting BEM solution\n')

%% Verification of the function - LAMBDA
% li=-0.2:0.001:0.4; fk=zeros(length(li),1);
% for k=1:length(li)
%         F=f(li(k),1);
%         fk(k)=F;
% end
% plot(li,fk)
% grid on

% % Verification of the function - PHI VALUES
% fi=-180:0.1:180; fk=zeros(length(fi),1); i=10; % Element
% thetaf=theta(i,1)*180/pi;
% 
% for k=1:length(fi)
%         [cl,cd]   = computeClCd(thetaf-fi(k), 0);         
%         F=f(0.05,fi(k));
%         fk(k)=F;
% end
% plot(fi,fk)
% grid on

%% NEWTON METHOD ITERATIONS
for j=1:length(Vc)
    fprintf('\nVc = %g m/s\n',Vc(j)) 
    thetaC_min = -10;     % angle de pas colectiu minim [deg] (aquell que fa que W>T)
    thetaC_max = 20;      % angle de pas colectiu maxim [deg] (aquell que fa que W<T)
    
    % Initialized here, no need to compute it every loop for finding theta
    Pi = 0;     % inicialitzacio de la potencia induida [W]
    P0 = 0;     % inicialitzacio de la potencia parasita [W] 
    Thrust = 0; % inicialitzaciï¿½ de l'empenta [N]
    
    for nt=1:nmaxITtheta
        
        % Array of colective command
        thetaC = 0.5*(thetaC_max+thetaC_min);   % calcul del punt mig de l'interval 
        theta_r= theta(:,j)*180/pi+thetaC;       % Vector de distribució de colectiu
        
        % Array de resultats intermitjos, a fi de veure la evolucio
        % iterativa
        li_verify=zeros(1,nmaxIT);
        fi_verify=zeros(1,nmaxIT);
        
        for i=1:nnodes
            
            %%% LAMBDA COMPUTATION - NEWTON METHOD
            %   Initial values
            li=0.05;
            fi=5;
            
            li_verify(1)=li;
            fi_verify(1)=fi;
            
            rf=r(i); sigf=sigma(i,j);
            
            for n=2:nmaxIT
                
                fi=atand((lc+li)/rf);
                [cl,cd]   = computeClCd(theta_r(i)-fi, 0);  
                [dcl,dcd] = compute_derClCd(theta_r(i)-fi, 0);
                  
                Ax=-Df(li,fi)\f(li,fi);
                li=li+Ax;        
                
                % Convergence Verification
                if max(abs([f(li,fi) fi-fi_verify(n-1)]))<prec_lambda % Result converged
                     lambda(i,j)=li;
                     phideg(i,j)=fi;
%                    fprintf('lambda_i at r(%g)=%g is equal to %g (n=%g)\n',i,r(i),li,n)
                    break
                elseif n==nmaxIT % Result did not converge
%                    fprintf('lambda_i at r(%g)=%g did not converge ---> mean approach\n',i,r(i))
%                    error('Amoavé')
                    li=mean(li_verify(50:100));
                    fi=mean(fi_verify(50:100));
                    lambda(i,j)=li;

%                 elseif i==10 % Review of the iteration process + results
%                    figure(3)
%                    subplot(1,2,1)
%                    plot(li_verify(1:n))
%                    subplot(1,2,2)
%                    plot(r,lambda(:,j))
%                    error('Amoavé')                    
                end 
                li_verify(n)=li;
                fi_verify(n)=fi;
            end
            [CL(i), CD(i)]=computeClCd(theta_r(i)-fi, 0); 
        end    
        % Verificació de correcte càlcul de velocitats induides
        plot(r,lambda(:,j))
    
    % Calcul del Thrust i les potencies
    K2 = CL.*cosd(phideg(:,j))-CD.*sind(phideg(:,j));  % coeficient utilitzat per partir l'expressiï¿½ del dT
    dT = -4*pi*rho*R^4*Omegadisseny^2.*lambda(:,j).*(lambda(:,j)+lambdac(j)).*r'*dr;
%    dFz = nb*0.5*rho*c(:,j)*vtip^2*R.*(r'.^2+(lambda(:,j)+lambdac(j)).^2).*K2*dr;   % Calucl diferencial de Thrust
    Thrust = sum(dT); 
    
    fprintf('IT %g: err=%g, T=%gN, W=%gN, thetaC=%gº\n',nt,abs(Thrust-W),Thrust,W,thetaC)
    %error('Please review the results so far for the Lambda solution')
    
        % Thrust convergence verification
        if abs(Thrust-W)<prec_thrust
            break
        elseif Thrust-W>=0
            thetaC_max = thetaC;        
        else
            thetaC_min = thetaC;                    
        end
    
        if nt==nmaxITtheta
            error('Theta no converge bro. A llorar a casa')
        end
    end
    
    % Guardem la solució
    thetaT(i,j) = theta(i,j)*(180/pi)+thetaC;    % angle theta [deg]     
    
end

fprintf('WIN')

%% Tornem a repetir-ho tot, pk podem
 for j=1:length(Vc)
    for i=1:nnodes

        phi_n(i,j) = atan(lambda(i,j)/r(i)); phideg(i,j) = phi_n(i,j)*(180/pi);            % angle de la velocitat incident [deg]
        alphadeg (i,j) = thetaT(i,j)-phideg(i,j);  alpha(i,j) = alphadeg(i,j)*(pi/180);                  % angle d'atac [deg]
        M = vtip*sqrt((lambda(i,j)+lambdac(j))^2+r(i)^2)/a; % numero de Mach a l'element de pala

        [cl(i,j),cd(i,j)] = computeClCd(alphadeg(i,j), M(i,j));               % funciï¿½ que dï¿½na Cl i Cd en funciï¿½ d'alpha i Mach

        K2 = cl(i,j)*cos(phi_n(i,j))-cd(i,j)*sin(phi_n(i,j));  % coeficient utilitzat per partir l'expressiï¿½ del dT

        dT(j) = nb*0.5*rho*c(i,j)*vtip^2*R*(r(i)^2+(lambda(i,j)+lambdac(j))^2)*K2*dr;   % Calucl diferencial de Thrust
        Thrust(j) = Thrust(j)+dT(j);            % Thrust de l'helicopter [N]
        
        dPi(j) = vtip*(lambda(i,j) + lambdac(j))*dT(j);   % Calcul diferencial de potencia induida
        Pi(j) = Pi(j)+dPi(j);                                 % Potencia induida
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
% title('BEM sense pÃ¨rdues');

