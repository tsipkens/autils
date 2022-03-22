
% MPFE_NI  Compute the mass-based filtration efficiency (or penetration).
%  Calculate via the numerical integration method.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-22

function [eta, s, G] = mpfe_ni(nup, ndown, di, prop, Gup, Gdown, szet)

Mup = pm.pm_ni(nup, di, prop);
Mdown = pm.pm_ni(ndown, di, prop);

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
    
    Nup = nansum(nup);
    Ndown = nansum(ndown);
    Eup = nansum(nup .* di .^ prop.zet);
    Edown = nansum(ndown .* di .^ prop.zet);
    
    J = P .* [-di .^ prop.zet ./ Eup; di .^ prop.zet ./ Edown; ...
        nansum(di .^ prop.zet .* log(di) .* (ndown ./ Edown - nup ./ Eup))];
    G = J' * G * J;  % LPU
    s = sqrt(diag(G));  % compute standard errors from covariance
else
    s = [];
    G = [];
end
%-------------------------------------------------------------------------%

end