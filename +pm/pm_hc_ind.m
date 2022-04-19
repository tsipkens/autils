
% PM_HC_IND  Compute the PM concentration using size-classified results with an independent N. 
%  Uncertainties in N should be appended as last entry in G.
%  
%  AUTHOR: Timothy Sipkens, 2022-04-19

function [M, s, G] = pm_hc_ind(Ni, di, prop, N, G)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop)  % if not given, assume water density and spheres
    prop = massmob.init('water');
end

di = di .* 1e9;

[dg, sg] = get_geo(Ni, di);

H = hc(dg, sg, prop.zet, prop.k);
M = H .* N;  % output mass

%-- UNCERTAINTIES --------------------------------------------------------%
if exist('G', 'var')
    % Detect if individual standard deviations are supplied. 
    % If so, convert to covariance matrix. 
    if any(size(G) == 1)
        G = diag(G .^ 2);
    end
    
    % Compute Jacobian. 
    Di = prop.zet .* (log(di) - log(dg)) + ...
        prop.zet ^ 2 ./ 2 .* ((log(di) - log(dg)) .^ 2 + log(sg) .^ 2);
    J = H .* Di;  % Jacobian
    
    J = [J; M ./ prop.rho100; ...
        prop.k * M * (prop.zet * log(sg) ^ 2 + log(dg)); ...
        H];  % supplement with mass-mob. parameters and N
    
    G = J' * G * J;  % LPU
    s = sqrt(G);  % compute standard errors from covariance
end
%-------------------------------------------------------------------------%

end

