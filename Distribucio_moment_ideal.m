%% Calcul distribucio del moment-ideal
dM=zeros(nnodes,length(Vc));
for i=1:nnodes
    for j=1:nVc
        for k=i:nnodes
            dM(i,j)=dM(i,j)+dL(k,j)*(r(k)-r(i))*R;
        end
    end
end
figure;
plot(r,dM);
xlabel('r')
ylabel('dM')
title('Distribucio de Moment - BEM ideal')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')