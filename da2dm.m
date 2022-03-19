
% DM2DA  Convert the aerodynamic diameter to a mobility diameter. 
%  
%  DM = da2dm(DA,RHO) converts the aerodynamic diameter, DA, to an 
%  mobility diameter, DM, using the particle density, RHO.
%  Expects diameters to be in m.
%  
%  DM = dm2da(DA,RHO,CHI) incorporates the dynamic shape factor, 
%  CHI. Default is CHI = 1. 
%  
%  DM = dm2da(DA,RHO,CHI,F_SIMPLE) adds a flag indicating whether 
%  to use the simple (F_SIMPLE = 1) or iterative (F_SIMPLE = 2) 
%  evaluation. The former is less accurate but faster. Default is
%  F_SIMPLE = 1.
%  
%  [DM,DVE] = da2dm(...) adds an output for the volume equivalent 
%  diameter of the particles. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2022-01-27

function [dm, dve] = da2dm(da, rho, chi, f_simple, varargin)

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
dve = da ./ sqrt(rho ./ rho0 ./ chi);  % volume equivalent diameter
dm = dve .* chi;  % mobility diameter

%-{
% Alternate, iterative method that is more precise.
opts = optimset('Display', 'off');
if ~f_simple
    % Solve for volume equivalent diameter.
    fun = @(dve) 1e9 .* (dve .* sqrt(rho ./ rho0 ./ chi .* ...
        Cc(dve, varargin{:}) ./ Cc(da, varargin{:})) - da);
    dve = fsolve(fun, dve, opts);
    
    % Solve for mobility diameter. 
    fun = @(dm) 1e9 .* (dm ./ chi .* ...
        Cc(dve, varargin{:}) ./ Cc(dm, varargin{:}) - dve);
    dm = fsolve(fun, dm, opts);
end
%}

end
