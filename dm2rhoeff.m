
% DM2RHOEFF  Get effective density from mobility diameter using mass-mob. relation.
%  
%  AUTHOR: Timothy Sipkens

function rho = dm2rhoeff(d, prop)

rho = rhoeff(d, dm2mp(d, prop));

end
