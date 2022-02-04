
% HC  Use Hatch-Choate to compute new distribution moments. 
%  
%  D = hc(MU, SG, Q) uses Hatch-Choate to compute new distribution 
%  moments having a GMD of MU and GSD of SG. Hatch-Choate applies the power
%  of Q. 
%  
%  D = hc(MU, SG, Q, A) adds a pre-factor A and switched to the integrated
%  variant. 
%  
%  AUTHOR: Timothy Sipkens, 2021-11-17

function d = hc(mu, sg, q, a)

if ~exist('a', 'var'); a = []; end

% Mean conversion.
d = mu .* exp(q .* (log(sg) .^ 2));

%-{
% Integrated variant if pre-factor, a, is specified. 
if ~isempty(a)
    d = a .* mu .^ q .* exp(1/2 .* q .^ 2 .* (log(sg) .^ 2));
end
%}

end

