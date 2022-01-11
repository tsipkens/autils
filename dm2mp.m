
% DM2MP  Calculate particle mass from a mobility diameter using mass-mobility relation.
%  
%  M = dm2mp(D,PROP) computes the particle mass from the mass-mobility
%  relation parameters in PROP (specifically PROP.Dm and prop.m0). 
% 
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function m = dm2mp(d, prop)

% Make sure 'm0' is a field of prop.
if ~isfield(prop, 'm0'); prop = prop_massmob(prop); end

% Use the mass-mobility relationship to get mobility diameter.
m = prop.m0 .* (d .* 1e9) .^ prop.Dm;

end

