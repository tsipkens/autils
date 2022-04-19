
% PM_HC_FIT  Calculate PM concentration using Hatch-Choate and lognormal fitting. 
%  
%  AUTHOR: Timothy Sipkens, 2022-04-19

function [M, s, G] = pm_hc_fit(Ni, di, prop, G)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop)  % if not given, assume water density and spheres
    prop = massmob.init('water');
end

di = di .* 1e9;

[dg, sg, Jg] = get_geo(Ni, di, 1);

H = hc(dg, sg, prop.zet, prop.k);
M = H .* nansum(Ni);  % output mass

%-- UNCERTAINTIES --------------------------------------------------------%
if exist('G', 'var')
    Gn = G(1:(end-2),1:(end-2));  % cov in Ni

    Gg = inv(Jg' * (Gn \ Jg));  % propagate through inverse procedure
    Gg = Gg(2:end, 2:end); % just dg and sg contributions
    
    % Jacobian and covariance for Hatch-Choate step.
    J = [H; M * prop.zet / dg; M * prop.zet ^ 2 * log(sg)];
    G = blkdiag(sum(sum(Gn)), Gg);
    G = J' * G * J;  % LPU
    
    s = sqrt(G);
end
%-------------------------------------------------------------------------%

end
