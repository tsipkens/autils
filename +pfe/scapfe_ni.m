
% SCAPFE_NI  Compute the scattering-based filtration efficiency (or penetration).
%  Calculate via the numerical integration method.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-26

function [eta, s, G] = scapfe_ni(nup, ndown, int)

P = nansum(ndown .* int) ./ nansum(nup .* int);

eta = 1 - P;

%-- UNCERTAINTIES --------------------------------------------------------%
% Not yet defined.
s = [];  G = [];
%-------------------------------------------------------------------------%

end
