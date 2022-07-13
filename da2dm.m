
% DM2DA  Convert the aerodynamic diameter to a mobility diameter. 
%  
%  DM = da2dm(DA,PROP) converts the aerodynamic diameter, DA, to an 
%  mobility diameter, DM, using the properties in PROP. 
%  
%  DM = dm2da(...,F_ITER) adds a flag indicating whether 
%  to use the simple (F_ITER = 0) or iterative (F_ITER = 1) 
%  evaluation. The former is less accurate but faster. Default is
%  F_SIMPLE = 1.
%  
%  [DM,DVE] = da2dm(...) adds an output for the volume equivalent 
%  diameter of the particles. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2022-01-27

function [dm, dve0] = da2dm(da, prop, f_iter, varargin)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end
%------------------------------------------%

rho0 = 1e3;  % density of water

% Compute simple volume-equivalent and aerodynamic diameters, 
% that is without iteration. 
opts = optimset('Display', 'off');
fun = @(dm) 1e9 .* ...
    (dm .* sqrt(dm2rhoeff(dm, prop) ./ rho0) - da);
dm = fsolve(fun, da, opts);

%-{
% Alternate, iterative method that is more precise.
if f_iter
    % Solve for mobility diameter. 
    fun = @(dm) 1e9 .* (dm .* ...
        sqrt(dm2rhoeff(dm, prop) ./ rho0 .* ...
        Cc(dm, varargin{:}) ./ Cc(da, varargin{:})) - da);
    dm = fsolve(fun, dm, opts);
end
%}

end
