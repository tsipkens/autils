
% DM_DA2MP  Uses mobility and aerodynamic diameter to compute mass.
%  
%  MP = dm_da2mp(DM, DA, VARARGIN) where DM is the mobility diameter [m], 
%  DA is the aerodynamic diameter [m], and the VARARGIN (optional) arguments  
%  are passed to the Cunningham slip correction calculation. 
%  
%  AUTHOR: Timothy Sipkens, 2024-01-12

function mp = dm_da2mp(dm, da, varargin)

% Compute simple volume-equivalent and aerodynamic diameters.
mp = da.^2 .* dm .* ...
    (pi*1e3/6 .* Cc(da, varargin{:}) ./ Cc(dm, varargin{:}));

end
