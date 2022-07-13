
% DVE2RHOEFF  Get effective density from volume-equivalent diameter using mass-mob. relation.
%  
%  AUTHOR: Timothy Sipkens

function rho = dve2rhoeff(dve, prop)

dm = dve2dm(dve, prop);
rho = rhoeff(dm, dm2mp(dm, prop));

end
