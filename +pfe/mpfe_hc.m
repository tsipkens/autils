
% MPFE_HC  Compute the mass-based filtration efficiency (or penetration).
%  Calculate via the Hatch-Choate method.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-22

function [eta, s, G] = mpfe_hc(nup, ndown, di, prop, Gup, Gdown, szet)

Mup = pm.pm_hc(nup, di, prop);
Mdown = pm.pm_hc(ndown, di, prop);

P = Mdown / Mup;
eta = 1 - P;

%-- UNCERTAINTIES --------------------------------------------------------%
% Detect if individual standard deviations are supplied. 
% If so, convert to covariance matrix. 
if exist('Gup', 'var')
    if ~exist('Gdown', 'var'); Gdown = []; end
    if ~exist('szet', 'var'); szet = 0; end
    
    if any(size(Gup) == 1)
        Gup = diag(Gup .^ 2);
    end
    if any(size(Gup) == 1)
        Gdown = diag(Gdown .^ 2);
    end
    G = blkdiag(Gup, Gdown, szet ^ 2);
    
    [dup, sup] = get_geo(nup, di);
    [ddown, sdown] = get_geo(ndown, di);

    Nup = nansum(nup);
    Ndown = nansum(ndown);
    Dup = 1 ./ Nup .* (1 + prop.zet .* (log(di) - log(dup)) .* ...
        prop.zet ^ 2 / 2 .* ((log(di) - log(dup)) .^ 2 - log(sup) ^ 2));
    Ddown = 1 ./ Ndown .* (1 + prop.zet .* (log(di) - log(ddown)) .* ...
        prop.zet ^ 2 / 2 .* ((log(di) - log(ddown)) .^ 2 - log(sdown) ^ 2));
    
    J = P .* [-Dup; Ddown; ...
        log(ddown) - log(dup) + (log(sup) ^ 2 - log(sdown) ^ 2)];
    G = J' * G * J;  % LPU
    s = sqrt(diag(G));  % compute standard errors from covariance
else
    s = [];
    G = [];
end
%-------------------------------------------------------------------------%

end