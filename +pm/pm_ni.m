
% PM_NI  Compute the PM concentration via numerical integration.
%  
%  AUTHOR: Timothy Sipkens, 04-2022

function [M, s, G, J] = pm_ni(Ni, di, prop, G)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop)  % if not given, assume water density and spheres
    prop.zet = 3;
    prop.rho100 = 1000;
end

di = di .* 1e9;  % mass-mobility relation define on diameters in nm

k = pi / 6 * prop.rho100 * (100^3) * 1e-27;
M = k .* nansum((di ./ 100) .^ prop.zet .* Ni);  % output mass
% M = prop.k .* nansum(di .^ prop.zet .* Ni, 1);  % pre-factor variant

%-- UNCERTAINTIES --------------------------------------------------------%
if exist('G', 'var')
    % Detect if individual standard deviations are supplied. 
    % If so, convert to covariance matrix. 
    if any(size(G) == 1)
        G = diag(G .^ 2);
    end
    
    % Compute Jacobian. 
    J = prop.k .* di .^ prop.zet;  % Jacobian
    if size(G, 1) == (length(di) + 2)  % then mass-mobility uncertainties are incl.
        J = [J; M ./ prop.rho100; ...
            k .* nansum((di ./ 100) .^ prop.zet .* Ni .* log(di ./ 100))];
        % J = [J; M ./ prop.k; prop.k .* ...
        %   nansum(di .^ prop.zet .* Ni .* log(di))];  % pre-factor variant
    end
    
    G = J' * G * J;  % LPU
    s = sqrt(G);  % compute standard errors from covariance
end
%-------------------------------------------------------------------------%

end
