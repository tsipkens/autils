
% MP2ZP  Calculate electric mobility from a vector of particle mass.
%  
%  B = mp2zp(M,Z,[],[],PROP) computes the mechanical mobility for the given particle
%  mass, M, and integer charge state, Z. Requires the PROP structure to 
%  have PROP.m0 and PROP.zet fields, which specify the mass-mobility
%  relation.
%  
%  B = mp2zp(M,Z,T,P,PROP) adds inputs explicitly stating the temperature
%  in Kelvin, T, and pressure in atm., P.
%  
%  [B,ZP] = mp2zp(...) add the electromobility, ZP, as an output. 
%  
%  [B,ZP,D] = mp2zp(...) add the mobility diameter implied by  the
%  mass-mobility relation.
% 
%  ------------------------------------------------------------------------
%  
%  NOTE:
%   Uses mass-mobility relationship to first convert to a mobility
%   diameter and then estimates the mobility using dm2zp.
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function [B, Zp, d] = mp2zp(m, z, T, P, prop)


%-- Parse inputs ---------------------------------------------------------%
if ~exist('T', 'var'); T = []; end
if ~exist('P', 'var'); P = []; end

% Make sure 'm0' is a field of prop.
if ~isfield(prop, 'm0'); prop = massmob.gen(prop); end
%-------------------------------------------------------------------------%


% Use the mass-mobility relationship to get mobility diameter.
d = 1e-9 .* (m ./ prop.m0) .^ (1 / prop.zet);

    
% Use mobility diameter to get particle 
% electrical and mechanical mobilities. 
if or(isempty(T), isempty(P))
    [B, Zp] = dm2zp(d, z);
else
    [B, Zp] = dm2zp(d, z, T, P);
end

end

