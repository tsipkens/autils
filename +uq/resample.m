
% RESAMPLE  Resample the point in Ni to generate new signals for UQ.
%  
%  AUTHOR: Timothy Sipkens, 2022-07-25

function [Nis] = resample(Ni, ns)

ni = size(Ni, 1);
n0 = size(Ni, 2);

idx2 = randi(n0, ni, ns);

Nis = zeros(ni, ns);
for ii=1:ni
    Nis(ii, :) = Ni(ii, idx2(ii, :));
end

end

