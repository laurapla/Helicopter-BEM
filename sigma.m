function s=sigma(alpha,r,lambda,theta)
s=(8*r*(lambda)^2)/((r^2+(lambda)^2)*(Cl(alpha)*cos(theta-(alpha*(pi/180)))-Cd(alpha)*sin(theta-(alpha*(pi/180)))));
end