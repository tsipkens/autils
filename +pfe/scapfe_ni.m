
% SCAPFE_NI  Compute the scattering-based filtration efficiency (or penetration).
%  Calculate via the numerical integration method.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-26

function eta = scapfe_ni(nup, ndown, di, l, the)

if ~exist('the', 'var'); the = []; end

if isempty(the)  % if no angle, use total scattering cross-section
    Qsca = mie.get_eff(l, di, 1.5442 + 0j);
    Csca = Qsca .* (pi .* di .^ 2 ./ 4);
    
else  % then use intensity at a specific angle
    % Actually intensity, Csca is a common surrogate var.
    Csca = mie.get_intensity(l, di, 1.5442 + 0j, [], the);
    
end

P = sum(ndown .* Csca) ./ sum(nup .* Csca);

eta = 1 - P;

end
