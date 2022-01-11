
% DM2DA  Convert the mobility diameter to an aerodynamic diameter. 
%  
%  DA = dm2da(DM,RHO) converts the mobility diameter, DM, to an 
%  aerodynamic diameter, DA, using the particle density, RHO. 
%  
%  DA = dm2da(DM,RHO,CHI) incorporates the dynamic shape factor, 
%  CHI. Default is CHI = 1. 
%  
%  DA = dm2da(DM,RHO,CHI,F_SIMPLE) adds a flag indicating whether 
%  to use the simple (F_SIMPLE = 1) or iterative (F_SIMPLE = 2) 
%  evaluation. The former is less accurate but faster. Default is
%  F_SIMPLE = 1.
%  
%  [DA,DVE] = dm2da(...) adds an output for the volume equivalent 
%  diameter of the particles. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function [da, dve] = dm2da(dm, rho, chi, f_simple, varargin)

%-- Parse inputs --------------------------%
if ~exist('rho', 'var'); rho = []; end
if isempty(rho); rho = 1e3; end

if ~exist('chi', 'var'); chi = []; end
if isempty(chi); chi = 1; end

if ~exist('f_simple', 'var'); f_simple = []; end
if isempty(f_simple); f_simple = 1; end
%-----------------------------------------%


rho0 = 1e3;  % density of water

% Compute simple volume-equivalent and aerodynamic diameters, 
% that is without iteration. 
dve = dm ./ chi;  % volume equivalent diameter
da = dve .* sqrt(rho ./ rho0 ./ chi);  % aerodynamic diameter

%-{
% Alternate, iterative method that is more precise.
opts = optimset('Display', 'off');
if ~f_simple
    fun = @(dve) 1e9 .* (dm ./ chi .* ...
        Cc(dve, varargin{:}) ./ Cc(dm, varargin{:}) - dve);
    dve = fsolve(fun, dve, opts);

    % Solve for aerodynamic diameter.
    fun = @(da) 1e9 .* (dve .* sqrt(rho ./ rho0 ./ chi .* ...
        Cc(dve, varargin{:}) ./ Cc(da, varargin{:})) - da);
    da = fsolve(fun, da, opts);
end
%}

end


