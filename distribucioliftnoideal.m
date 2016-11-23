figure;
plot (r,dLPrandtl0,r,dLPrandtl25,r,dLPrandtl5,r,dLPrandtl75,r,dLPrandtl10,r,dLPrandtl125)
xlabel('r')
ylabel('dL')
title('Distribucio de Lift - BEM+perdues+compressibilitat')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')

dM0=dLPrandtl0.*r*R;
dM25=dLPrandtl25.*r*R;
dM5=dLPrandtl5.*r*R;
dM75=dLPrandtl75.*r*R;
dM10=dLPrandtl10.*r*R;
dM125=dLPrandtl125.*r*R;

figure;
plot (r,dM0,r,dM25,r,dM5,r,dM75,r,dM10,r,dM125)
xlabel('r')
ylabel('dM')
title('Distribucio de Moments - BEM+perdues+compressibilitat')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')