
% MC  Monte Carlo method.
%  
%  AUTHOR: Timothy Sipkens, 2022-02-10

function [y, s, smpl, f_smpl] = mc(x, Gx, f, n, seed)

if ~exist('seed', 'var'); seed = []; end
if isempty(seed); seed = randi(1); end
rng(seed);  % reset number generator with seed

% If vector of variances supplied.
if any(size(Gx) == 1)
    Gx = diag(Gx .^ 2);
end


smpl = mvnrnd(x, Gx, n)';

f_smpl = f(smpl);
y = mean(f_smpl, 2);
s = std(f_smpl, [], 2);

end
