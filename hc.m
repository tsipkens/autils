
% HC  Use Hatch-Choate to compute new distribution moments. 
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

