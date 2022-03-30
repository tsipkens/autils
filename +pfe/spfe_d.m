
% SPFE_D  Compute the PFE at a specific size (e.g., PFE at 100 nm).
%  
%  NOTE: A single covariance matrix can be provided, allowing for
%  correlation between upstream and downstream locations. 
%  
%  AUTHOR: Timothy Sipkens, 2022-03-22

function [eta_d, s_d] = spfe_d(nup, ndown, di, dq, Gup, Gdown)

% First, get size-resolved PFEs.
if exist('Gup', 'var')
    [eta, s] = pfe.spfe(nup, ndown, Gup, Gdown);
else
    eta = pfe.spfe(nup, ndown);
end

% Then interpolate, using loop for each input set.
eta_d = zeros([1, size(nup, 2)]);
if exist('Gup', 'var')
    s_d = zeros([1, size(nup, 2)]);
end
for ii=1:size(nup, 2)
    eta_d(:,ii) = interp1(di(:,ii), eta(:,ii), dq);  % interpolate for value at a specific size
    
    if exist('Gup', 'var')
        s_d(:,ii) = interp1(di(:,ii), s(:,ii), dq);  % blind interpolation of uncertainties (not the best)
    end
end

end
