
% MC  Monte Carlo method.
%  F_NEG determines if non-negativity is enforced for the samples.
%  
%  AUTHOR: Timothy Sipkens, 2022-02-10

function [y, s, smpl, f_smpl] = mc(x, Gx, f, n, seed, f_neg)

if ~exist('f_neg', 'var'); f_neg = 1; end
if isempty(f_neg); f_neg = 1; end

if ~exist('seed', 'var'); seed = []; end
if isempty(seed); seed = randi(1e5); end
rng(seed);  % reset number generator with seed

% If vector of variances supplied.
if any(size(Gx) == 1)
    Gx = diag(Gx .^ 2);
end


smpl = mvnrnd(x, Gx, n)';
if f_neg; smpl(:, any(smpl < 0)) = []; end

f_smpl = f(smpl);
y = mean(f_smpl, 2);
s = std(f_smpl, [], 2);

end
