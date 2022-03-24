
% MP2DM  Calculate particle mass from a mobility diameter using mass-mobility relation.
%  
%  D = mp2dm(M,PROP) computes the mobility diameter from the mass-mobility
%  relation parameters in PROP (specifically PROP.Dm and prop.m0). 
% 
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function d = mp2dm(m, prop)

% Make sure 'm0' is a field of prop.
if ~isfield(prop, 'm0'); prop = massmob.init(prop); end

% Use the mass-mobility relationship to get mobility diameter.
d = 1e-9 .* (m ./ prop.m0) .^ (1 / prop.zet);

end

