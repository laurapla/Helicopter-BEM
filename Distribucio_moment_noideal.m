%% Distribucio_moments_noideal
for i=1:nnodes
    for j=1:length(Vc)
        dM(i,j)=dL(i,j)*r(i)*R;
    end
end

figure;
plot(r,dM);
xlabel('r')
ylabel('dM')
title('Distribucio de Moments - BEM sense perdues')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')