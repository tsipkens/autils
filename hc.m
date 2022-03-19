
% HC  Use Hatch-Choate to compute new distribution moments. 
%  
%  D = hc(MU, SG, Q) uses Hatch-Choate to compute new distribution 
%  moments for a distribution having a GMD of MU and GSD of SG, using a
%  power of Q (e.g., the mass-mobility exponent).
%  
%  D = hc(MU, SG, Q, A) adds a pre-factor A and switches to the integrated
%  variant, which is useful for computing total PM or number concentration.
%  
%  AUTHOR: Timothy Sipkens, 2021-11-17

function d = hc(mu, sg, q, a)

if ~exist('a', 'var'); a = []; end

if isempty(a)  % Transform the mean.
    d = mu .* exp(q .* (log(sg) .^ 2));

else % Integrated variant if pre-factor, a, is specified. 
    d = a .* mu .^ q .* exp(1/2 .* q .^ 2 .* (log(sg) .^ 2));
end

end

