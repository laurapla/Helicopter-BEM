%% Distribucio_lift_noideal
for i=1:nnodes
    for j=1:length(Vc)
        dL(i,j)=0.5*rho*c(i,j)*sqrt((Omegadisseny*R*r(i))^2+(Vc(j)+lambda(i,j)*Omegadisseny*R)^2)*dr*R;
    end
end

figure;
plot(r,dL);
xlabel('r')
ylabel('dL')
title('Distribucio de Lift - BEM sense perdues')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')