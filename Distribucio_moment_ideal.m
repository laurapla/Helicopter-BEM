%% Calcul distribucio del moment-ideal
for i=1:nnodes
    for j=1:nVc
        dM(i,j)=dL(i,j)*r(i)*R;
    end
end
figure;
plot(r,dM);
xlabel('r')
ylabel('dM')
title('Distribucio de Moment - BEM ideal')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')