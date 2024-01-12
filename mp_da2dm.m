
% MP_DA2DM  Uses particle mass and aerodynamic diameter to compute mobility diameter.
%  
%  DM = mp_da2dm(MP, DA, F_ITER, VARARGIN) where MP is the particle mass 
%  [kg], DA is the aerodynamic diameter [m], F_ITER (optional) flags 
%  whether to apply an iterative method to correct for slip, and the 
%  VARARGIN (optional) arguments are passed to the Cunningham slip  
%  correction calculation. 
%  
%  AUTHOR: Timothy Sipkens, 2024-01-12

function dm = mp_da2dm(mp, da, f_iter, varargin)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end
%-----------------------------------------%

% Compute simple volume-equivalent and aerodynamic diameters.
dm = da.^(-2) .* mp .* (6/(pi*1e3));

%-{
% Alternate, iterative method that is more precise.
opts = optimset('Display', 'off');
if f_iter
    % Solve for aerodynamic diameter.
    fun = @(dm1) 1e9 .* ...
        (dm .* (Cc(dm1, varargin{:}) ./ Cc(da, varargin{:})) - dm1);
    dm = fsolve(fun, dm, opts);
end
%}

end
