
% EQUIV_DIST  An adaptive function to determine equivalent lognormal distribution properties.
%  
%  By default, takes mobility diameter parameters as inputs and converts
%  to aerodynamic diameter, volume-equivalent diameter, and mass.
%  
%  AUTHOR: Timothy Sipkens, 2022-04-14

function tb = equiv_dist(cmd, sg, prop, type)

if ~exist('type', 'var'); type = ''; end

% type = 'dm';
f_aero = 0;

if strcmp(type, 'mp')
    sg = exp( ...
        (log(mp2dm(exp(log(cmd) .* 1.05), prop)) - ...
         log(mp2dm(exp(log(cmd) .* 0.95), prop))) ./ ...
        (0.1 .* log(cmd)) .* ...
        log(sg));  % mass
    
    cmd = mp2dm(cmd, prop);
end

rho = dm2rhoeff(cmd, prop);
chi = prop.chi;

% Transform CMDs. 
[cmad, cmvd] = dm2da(cmd, rho, chi, f_aero);
mg = dm2mp(cmd, prop);

% Transform GSDs.
[dau, vdu] = dm2da(exp(log(cmd) .* 1.05), rho, chi, f_aero);
[dad, vdd] = dm2da(exp(log(cmd) .* 0.95), rho, chi, f_aero);

sga = exp( ...
    (log(dau) - ...
     log(dad)) ./ ...
    (0.1 .* log(cmd)) .* ...
    log(sg)); % aerodynamic diameter

sgv = exp((log(vdu) - log(vdd)) ./ ...
    (0.1 .* log(cmd)) .* ...
    log(sg)); % aerodynamic diameter

sgm = exp( ...
    (log(dm2mp(exp(log(cmd) .* 1.05), prop)) - ...
     log(dm2mp(exp(log(cmd) .* 0.95), prop))) ./ ...
    (0.1 .* log(cmd)) .* ...
    log(sg));  % mass

% Get table of values.
s.Mobility = [cmd .* 1e9; sg];
s.Aerodynamic = [cmad .* 1e9; sga];
s.Volume = [cmvd .* 1e9; sgv];
s.Mass = [mg .* 1e18; sgm];

tb = struct2table(s);
tb.Properties.RowNames = {'GMD', 'GSD'};

end
