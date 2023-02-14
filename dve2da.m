
% DVE2DA  Convert the mobility diameter to an aerodynamic diameter. 
%  
%  DA = dve2da(DVE,PROP) performs the conversion using the  
%  mass-mobility properties and material density in PROP.
%  
%  DA = dve2da(DVE,PROP,F_ITER) adds a flag indicating whether 
%  to use the simple (F_ITER = 0) or iterative (F_ITER = 1) 
%  evaluation. The former is less accurate but faster. 
%  Default is F_ITER = 1.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function da = dve2da(dve, prop, f_iter, varargin)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end  % default is iterative method
%-----------------------------------------%

rho0 = 1e3;  % density of water

% Compute simple volume-equivalent and aerodynamic diameters. 
% that is without iteration. 
chi = dve2chi(dve, prop, f_iter);
da = dve .* sqrt(prop.rhom ./ rho0 ./ chi);

%-{
% Alternate, iterative method that is more precise.
opts = optimset('Display', 'off');
if f_iter
    % Solve for aerodynamic diameter.
    fun = @(da) 1e9 .* (dve .* sqrt(prop.rhom ./ rho0 ./ chi .* ...
        Cc(dve, varargin{:}) ./ Cc(da, varargin{:})) - da);
    da = fsolve(fun, da, opts);
end
%}

end
