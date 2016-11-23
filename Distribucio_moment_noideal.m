%% Distribucio_moments_noideal
dM=zeros(nnodes,length(Vc));
for i=1:nnodes
    for j=1:length(Vc)
        for k=i:nnodes
            dM(i,j)=dM(i,j)+dL(k,j)*(r(k)-r(i))*R;
        end
    end
end

figure;
plot(r,dM);
xlabel('r')
ylabel('dM')
title('Distribucio de Moments - BEM sense perdues')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')