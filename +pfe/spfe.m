
% SPFE  Compute the size-resolved filtration efficiency (or penetration).
%  
%  NOTE: A single covariance matrix can be provided, allowing for
%  correlation between upstream and downstream locations. 
%  
%  AUTHOR: Timothy Sipkens, 2022-03-22

function [eta, s, G] = spfe(nup, ndown, Gup, Gdown, f_sampl)

if ~exist('f_sampl', 'var'); f_sampl = []; end
if isempty(f_sampl); f_sampl = 0; end

P = ndown ./ nup;
eta = 1 - P;

%-- UNCERTAINTIES --------------------------------------------------------%
% Detect if individual standard deviations are supplied. 
% If so, convert to covariance matrix. 
if exist('Gup', 'var')
    if ~exist('Gdown', 'var'); Gdown = []; end
    
    % Requires standard errors, not covariances. 
    % Then, Monte Carlo sampling is used to compute intervals. 
    % This is more important for SPFEs, where lower counts may be 
    % encountered at the edges of the size distribution.
    if f_sampl
        ns = 400;  % number of samples
        
        % Sample quantities, with repeats in third dimension.
        nsup = normrnd(repmat(nup, [1,1,ns]), repmat(Gup, [1,1,ns]));
        nsdown = normrnd(repmat(ndown, [1,1,ns]), repmat(Gdown, [1,1,ns]));
        
        s = std(nsdown ./ nsup, [], 3);  % standard deviation on thrid dim.
        G = [];  % no covariance in this treatment
        
    else
        % If a single set of size-resolved measurements, allow for cov.
        if any(size(nup) == 1)
            if any(size(Gup) == 1)
                Gup = diag(Gup .^ 2);
            end
            if any(size(Gdown) == 1)
                Gdown = diag(Gdown .^ 2);
            end
            G = blkdiag(Gup, Gdown);

            J = [diag(-P ./ nup); diag(P ./ ndown)];  % Jacobian

            G = J' * G * J;  % LPU
            s = sqrt(diag(G));  % compute standard errors from covariance


        % If a set of multiple size-resolved measurements. 
        % Only considers standard error (i.e., no cov), 
        % which simplifies code and input variables.
        else
            s = sqrt((P ./ nup .* Gup) .^ 2 + (P ./ ndown .* Gdown) .^ 2);
            G = [];  % don't output this
        end
    end
else
    s = [];
    G = [];
end
%-------------------------------------------------------------------------%

end
