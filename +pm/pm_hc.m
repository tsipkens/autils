

function [M, s, G, J] = pm_hc(Ni, di, prop, G)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop)  % if not given, assume water density and spheres
    prop = massmob.init('water');
end

di = di .* 1e9;

[dg, sg] = get_geo(Ni, di);

H = hc(dg, sg, prop.zet, prop.k);
M = H .* nansum(Ni);  % output mass

%-- UNCERTAINTIES --------------------------------------------------------%
if exist('G', 'var')
    % Detect if individual standard deviations are supplied. 
    % If so, convert to covariance matrix. 
    if any(size(G) == 1)
        G = diag(G .^ 2);
    end
    
    % Compute Jacobian. 
    Di = 1 + prop.zet .* (log(di) - log(dg)) + ...
        prop.zet ^ 2 ./ 2 .* ((log(di) - log(dg)) .^ 2 + log(sg) .^ 2);
    J = H .* Di;  % Jacobian

    %{
    if size(G, 1) == (length(di) + 2)  % then mass-mobility uncertainties are incl.
        J = [J; M ./ prop.rho100;  ...
            prop.k * M * (prop.zet * log(sg) ^ 2 + log(dg))];
    end
    %}
    if size(G, 1) == (length(di) + 2)  % then mass-mobility uncertainties are incl.
        J = [J; M ./ prop.rho100;  ...
            M * (prop.zet * log(sg) ^ 2 + log(dg ./ 100))];
    end
    
    G = J' * G * J;  % LPU
    s = sqrt(G);  % compute standard errors from covariance
end
%-------------------------------------------------------------------------%

end
