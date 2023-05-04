
% RGD2DP  Convert a radius of gyration to primary particle using dimensions/prefactors.
%  
%  Takes radius of gyration as an input and uses the following dimensions
%  and prefactors in calculation to bypass need for explicit Np: 
%  kf, Df, dp100, Dtem, k_alpha, D_alpha. 
%  
%  Derivation involves using relationship to mobility diameter, but
%  substituting this quantity in two equations to avoid its direct
%  calculation. 
%  
%  AUTHOR: Timothy Sipkens, 2023-04-20

function [dp, Np] = rgd2dp(Rg, prop)

dp = (prop.kf ./ prop.k_alpha .* (2 .* Rg) .^ prop.Df .* ...
    (100e-9 ./ prop.dp100 .^ (1/prop.Dtem)) .^ (-2 * prop.D_alpha)) .^ ...
    (1 ./ (prop.Df + 2 * prop.D_alpha * (1/prop.Dtem - 1)));

Np = prop.kf .* (2 .* Rg ./ dp) .^ prop.Df;

end
