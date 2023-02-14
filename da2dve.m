
% DA2DVE  Convert the aerodynamic diameter to a volume-equivalent diameter. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2022-07-13

function dve = da2dve(da, prop, f_iter)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end
%------------------------------------------%

rho0 = 1e3;

if (f_iter)
    fun = @(dve) (dve .* sqrt(prop.rhom / rho0 ./ dve2chi(dve, prop, f_iter) * Cc(dve) / Cc(da)) - da) .* 1e9;
else
    fun = @(dve) (da ./ sqrt(prop.rhom / rho0 ./ dve2chi(dve, prop, f_iter)) - dve) .* 1e9;
end

opts = optimset('Display', 'off');
dve = fsolve(fun, da, opts);

end
