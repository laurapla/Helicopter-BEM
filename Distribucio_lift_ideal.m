%% Calcul de la distribucio de Lift-ideal

nelem=1000;
nnodes=nelem+1;

%% BEM ideal (MTH)


% Discretització
r = linspace(0.1,1,nnodes);
dr = r(end)/nelem;

% Càlcul de angles i solidesa
% Cada columna de la matriu conte la variable calculada per a una velocitat de climbing
sigmaideal = zeros(nnodes,nVc);
cideal = zeros(nnodes,nVc);
nb = 6; % calculat més endavant amb la BEM
Omegadisseny = vtip/R; %[rad/s]

for j=1:nVc  % Calcula la solides i angles per tots els valors de Vc
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
for i=1:nnodes
    for j=1:nVc
        dL(i,j)=0.5*rho*Cl*cideal(i)*((Omegadisseny*r(i)*R)^2+(Vc(j)+vi)^2)*dr*R;
    end
end

figure;
plot(r,dL);
xlabel('r')
ylabel('dL')
title('Distribucio de Lift - BEM ideal')

