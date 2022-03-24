
% RHO_EFF  Computes the effective density from mass and mobility diameter.
%  
%  Expects mobility diameter in m and mass in kg.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-18

function rho = rho_eff(dm, m)

rho = 6 .* m ./ (pi .* dm .^ 3);

end
