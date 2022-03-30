
% SCAPFE_NI  Compute the scattering-based filtration efficiency (or penetration).
%  Calculate via the numerical integration method.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-26

function [eta, s, G] = scapfe_ni(nup, ndown, di, prop, Gup, Gdown, szet)

P = sum(ndown .* di .^ 6) ./ sum(nup .* di .^ 6);

eta = 1 - P;

end

