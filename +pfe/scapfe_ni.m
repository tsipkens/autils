
% SCAPFE_NI  Compute the scattering-based filtration efficiency (or penetration).
%  Calculate via the numerical integration method.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-26

function eta = scapfe_ni(nup, ndown, di, l)

Qsca = mie.get_mie(l, di, 1.5442 + 0j);
Csca = Qsca .* (pi .* di .^ 2 ./ 4);

P = sum(ndown .* Csca) ./ sum(nup .* Csca);

eta = 1 - P;

end

