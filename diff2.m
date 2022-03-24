
% DIFF2  Computes numerical derivative in log-space given two vectors.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-24

function [d, dd1, dd2] = diff2(d1, d2)

dd1 = log(d1(2:end)) - log(d1(1:end-1));
dd1 = [dd1, dd1(end)];

dd2 = log(d2(2:end)) - log(d2(1:end-1));
dd2 = [dd2, dd2(end)];

d = dd1 ./ dd2;

end
