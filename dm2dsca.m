
% DM2DSCA  Convert mobility diameter to scattering-equivalent diameter.
%  
%  AUTHOR: Timothy Sipkens, 2021-11-17

function [dsca] = dm2dsca(d, m, k)

F = abs((m .^ 2 - 1) ./ (m .^ 2 + 2)) .^ 2;

sig = (8*pi/3) .* k .^ 4 .* (d/2) .^ 6 .* F;  % scattering cross section
dsca = (sig ./ (pi/6)) .^ (1/3);

end

