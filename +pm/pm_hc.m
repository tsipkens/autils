

function [M, G] = pm_hc(Ni, di, prop, G)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop)  % if not given, assume water density and spheres
    prop.zet = 3;
    prop.k = 1000 * pi / 6;
end

di = di .* 1e9;

[dg, sg] = get_geo(Ni, di);

k = pi / 6 * prop.rho100 * (100^3) * 1e-27;
H = hc(dg / 100, sg, prop.zet, k);
M = H * nansum(Ni);  % output mass

%-- UNCERTAINTIES --------------------------------------------------------%
% Detect if individual standard deviations are supplied. 
% If so, convert to covariance matrix. 
if any(size(G) == 1)
    G = diag(G .^ 2);
end

% Compute Jacobian. 
Di = 1 + prop.zet .* (log(di) - log(dg)) + ...
    prop.zet ^ 2 .* ((log(di) - log(dg)) .^ 2 + log(sg) .^ 2);
J = H .* Di;  % Jacobian
if size(G, 1) == (length(di) + 2)  % then mass-mobility uncertainties are incl.
    J = [J; M ./ prop.rho100; k .* ...
        M * (prop.zet * log(sg) ^ 2 + log(dg))];
end

G = J' * G * J;  % LPU
%-------------------------------------------------------------------------%

end
