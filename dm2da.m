
% DM2DA  Convert the mobility diameter to an aerodynamic diameter. 
%  
%  DA = dm2da(DM,PROP) converts the mobility diameter, DM, to an 
%  aerodynamic diameter, DA, using the properties in PROP.
%  Expects diameters to be in m.
%  
%  DA = dm2da(...,F_ITER) adds a flag indicating whether 
%  to use the simple (F_ITER = 0) or iterative (F_ITER = 2) 
%  evaluation. The former is less accurate but faster. 
%  Default is F_ITER = 1. 
%  
%  NOTE: Does not require PROP.rhom, which is avoided by using 
%  an "representative" volume-equivalent diameter that ignores the 
%  shape factor.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function [da, dve0] = dm2da(dm, prop, f_iter, varargin)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end
%-----------------------------------------%

rho0 = 1e3;  % density of water

% Compute simple volume-equivalent and aerodynamic diameters, 
% that is without iteration. 
rho = dm2rhoeff(dm, prop);  % effective density
da = dm .* sqrt(rho ./ rho0);  % aerodynamic diameter (simple)

%-{
% Alternate, iterative method that is more precise.
opts = optimset('Display', 'off');
if f_iter
    % Solve for aerodynamic diameter.
    fun = @(da) 1e9 .* (dm .* sqrt(rho ./ rho0 .* ...
        Cc(dm, varargin{:}) ./ Cc(da, varargin{:})) - da);
    da = fsolve(fun, da, opts);
end
%}

end
