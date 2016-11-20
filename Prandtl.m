function F = Prandtl(phi,nb,r,rroot)

% This function computes the F factor (Prandtl)
% - phi [rad]
% - nb: number of blades
% - r: radial position - adimensional (x/R)
% - rroot: radial postion of the first point of the blade (not hub) -
%   adimensional (x/R)

ftip = nb/2*((1-r)/(r*sin(phi)));
froot = nb/2*((r-rroot)/(r*sin(phi)));
F = 4/pi^2*acos(exp(-ftip))*acos(exp(-froot));

end