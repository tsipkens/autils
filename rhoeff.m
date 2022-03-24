
% RHOEFF  Computes the effective density from mass and mobility diameter.
%  
%  RHO = rhoeff(D,M) computes the effective density using the mobility
%  diameter, D, and particle mass, M. Expects the mobility diameter, d, 
%  to be in m and the mass to be in kg.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-18

function rho = rhoeff(d, m)

rho = 6 .* m ./ (pi .* d .^ 3);

end
