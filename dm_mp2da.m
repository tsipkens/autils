
% DM_MP2DA  Uses mobility diameter and mass to compute aerodynamic diameter.
%  
%  DA = dm_mp2da(DM, MP, F_ITER, VARARGIN) where DM is the mobility  
%  diameter [m], MP is the particle mass [fg], F_ITER (optional)  
%  flags whether to apply an iterative method to correct for slip, and the 
%  VARARGIN (optional) arguments are passed to the Cunningham slip  
%  correction calculation. 
%  
%  AUTHOR: Timothy Sipkens, 2024-01-12

function da = dm_mp2da(dm, mp, f_iter, varargin)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end
%-----------------------------------------%

% Compute simple volume-equivalent and aerodynamic diameters.
da = dm.^(-1/2) .* mp.^(1/2) .* sqrt(6/(pi*1e3));

%-{
% Alternate, iterative method that is more precise.
opts = optimset('Display', 'off');
if f_iter
    % Solve for aerodynamic diameter.
    fun = @(da1) 1e9 .* ...
        (da .* sqrt(Cc(dm, varargin{:}) ./ Cc(da1, varargin{:})) - da1);
    da = fsolve(fun, da, opts);
end
%}

end
