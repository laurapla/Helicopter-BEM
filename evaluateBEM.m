function [F] = evaluateBEM(L,Lc,r,theta,sigma)
 
    phi = atand((L+Lc)/r);               % angle de la velocitat incident [deg]
    alpha = theta-phi;                              % angle d'atac [deg]
    [cl,cd] = computeClCd(alpha, 0);                       % coeficients aerodinamics de l'element de pala

    % Equacio de segon grau en lambda induida per l'element i
    K = (cl*cosd(phi)-cd*sind(phi))*sigma;

    F=8*r*(Lc+L)*L-K*(r^2+(L^2+Lc^2));

end