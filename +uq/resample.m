
% RESAMPLE  Resample the point in Ni to generate new signals for UQ.
%  
%  AUTHOR: Timothy Sipkens, 2022-07-25

function [Nis] = resample(Ni, ns)

ni = size(Ni, 1);
n0 = size(Ni, 2);

idx2 = randi(n0, ni, ns);
idx1 = (1:ni)' * ones(1, ns);

idx = sub2ind(size(Ni), idx1, idx2);

Nis = Ni(idx);

end

