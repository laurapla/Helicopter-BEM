syms Vi
vi=zeros(1,6);
for i=1:nVc
    G = double(solve(W==2*rho*A*(Vc(i)+Vi)*Vi,Vi));
    vi(1,i)=G(2);
    Omega(i)=(sqrt(vtip^2-(Vc(i)+vi(1,i))^2))/R;
end
lambdac = Vc / vtip;

for i=1:length(r)
    for j=1:nVc
        lambdai(i,j)=(sqrt(vtip^2-(Omega(j)*R)^2)-Vc(j))/vtip;
    end
end
figure;
plot(r,lambdai);
xlabel('r')
ylabel('\lambda_{i}')
title('Velocitat induïda - BEM ideal')
legend('Vc=0','Vc=2,5m/s','Vc=5m/s','Vc=7,5m/s','Vc=10m/s','Vc=12,5m/s')
