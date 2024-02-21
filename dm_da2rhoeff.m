
% DM_DA2RHOEFF  Uses mobility and aerodynamic diameter to compute eff. dens.
%  
%  RHO = dm_da2rhoeff(DM, DA, ...) uses mobility, DM, and aerodynamic, DA,
%  diameters to compute effective density. Additional arguments are passed
%  to the slip correction. 
%  
%  AUTHOR: Timothy Sipkens, 2024-02-21

function rho = dm_da2rhoeff(dm, da, varargin)

rho = 1000 .* (da ./ dm) .^2 .* Cc(da, varargin{:}) ./ Cc(dm, varargin{:});

end
