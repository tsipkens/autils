
% ZP2DM  Calculate mobility diameter from a vector of mobilities.
% 
%  B = zp2dm(DM,Z) computes the mechanical mobility for the given particle
%  mobility diameter, D, and integer charge state, Z. 
%  DM is to be given in m. 
%  
%  B = zp2dm(DM,Z,T,P) adds inputs explicitly stating the temperature
%  in Kelvin, T, and pressure in atm., P.
%  
%  [B,ZP] = mp2zp(...) add the electromobility, ZP, as an output. 
%  
%  ------------------------------------------------------------------------
%  
%  NOTES:
%   1 Temperature (T) and pressure (p) must both be specified before 
%     they are used. Otherwise, they are ignored and the code from Buckley
%     et al. (2017) is used.
%   2 Some of the code is adapted from Buckley et al. (2017) and Olfert 
%     laboratory.
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function [dm, B] = zp2dm(Zp, z, T, p)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('z','var'); z = []; end
if isempty(z); z = 1; end  % if integer charge state not specified use z = 1

%-- Perform calculation --------------------------------------------------%
e = 1.6022e-19; % define electron charge [C]
B = Zp ./ (e .* z); % mechanical mobility

if nargin <= 3  % if P and T are not specified, use Buckley/Davies
    mu = 1.82e-5;  % gas viscosity [Pa*s]
    Cc_fun = @(dm) Cc(dm);
    
else % If P and T are Olfert laboratory / Kim et al.
    S = 110.4; % temperature [K]
    T_0 = 296.15; % reference temperature [K]
    vis_23 = 1.83245e-5; % reference viscosity [kg/(m*s)]
    mu = vis_23 * ((T / T_0) ^ 1.5) * ((T_0 + S)/(T + S)); % gas viscosity
        % Kim et al. (2005), ISO 15900, Eqn 3

    Cc_fun = @(dm) Cc(dm, T, p);
end

dm0 = 1 ./ (3 * pi * mu .* B); % initialize mechanical mobility
dm0 = dm0 .* Cc_fun(dm0);

opts = optimset('Display', 'off');
dm = fminsearch(@(dm) ...
    1e9 .* norm(B - Cc_fun(abs(dm)) ./ (3 * pi * mu .* abs(dm))), dm0, opts);
%-------------------------------------------------------------------------%

end


