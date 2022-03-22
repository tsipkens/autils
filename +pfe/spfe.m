
% SPFE  Compute the size-resolved filtration efficiency (or penetration).
%  
%  NOTE: A single covariance matrix can be provided, allowing for
%  correlation between upstream and downstream locations. 
%  
%  AUTHOR: Timothy Sipkens, 2022-03-22

function [eta, s, G] = spfe(nup, ndown, Gup, Gdown)

P = ndown ./ nup;
eta = 1 - P;

%-- UNCERTAINTIES --------------------------------------------------------%
% Detect if individual standard deviations are supplied. 
% If so, convert to covariance matrix. 
if exist('Gup', 'var')
    if ~exist('Gdown', 'var'); Gdown = []; end
    
    if any(size(Gup) == 1)
        Gup = diag(Gup .^ 2);
    end
    if any(size(Gup) == 1)
        Gdown = diag(Gdown .^ 2);
    end
    G = blkdiag(Gup, Gdown);
    
    J = [diag(-P ./ nup); diag(P ./ ndown)];  % Jacobian

    G = J' * G * J;  % LPU
    s = sqrt(diag(G));  % compute standard errors from covariance
else
    s = [];
    G = [];
end
%-------------------------------------------------------------------------%

end
