
% SDM2SDA  Calculate the GSD for aerodynamic distribution from mobility distribution GSD.
%  
%  SM = sdm2sda(SD,PROP) converts the GSD using the mass-mobility
%  relation parameters in PROP (specifically PROP.zet). 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2022-07-13

function sa = sdm2sda(sd, d, prop, f_iter)

if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 0; end

% Perturb to estimate new GSD.
d1 = exp(log(d) .* 1.01);
d2 = exp(log(d) .* 0.99);
sa = exp( ...
    (log(dm2da(d1, prop, f_iter)) - ...
     log(dm2da(d2, prop, f_iter))) ./ ...
    (0.02 .* log(d)) .* ...
    log(sd));  % aerodynamic diameter, count distribution GSD

end
