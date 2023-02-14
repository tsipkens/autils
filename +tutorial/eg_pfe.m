
% Example of computing uncertainties for PFE.

clear;
close all;


%== GENERATE NOISY SIGNALS ===============================================%
%   Can be replaced with experimental measurements.
rng(10);  % control random number generator
eta = 0.999;  % pfe < CHANGE THIS
P = 1 - eta;  % penetration
N0 = 1e6;  % total number concentration upstream
n = 100;  % number of samples

% Add Gaussian noise.
s1 = 0.4 * N0;  % < CHANGE THIS
s2 = 0.4 * P * N0;  % < CHANGE THIS
N1 = normrnd(N0, s1, [1, n]);  % upstream "CPC" counts
N2 = normrnd(P .* N0, s2, [1, n]);  % downstream "CPC" counts
Ni = [N1; N2];  % compile inputs
%=========================================================================%


n = size(Ni, 2);  % compute number of repeats (if experimental signals replace above)

%== STEP 0: Covariances and expected values ==%
%   These are common to the difference methods below.
Ni_bar = mean(Ni, 2);  % get mean of measurements
G = cov(Ni');  % get covariance


%== LPU ==================================================================%
%-- STEP 0: Covariances and expected values --%
P_lpu = Ni_bar(2) ./ Ni_bar(1); % expected penetration

%-- STEP 2: Jacobian --%
J = P_lpu .* [1 ./ Ni_bar(1); 1 ./ Ni_bar(2)];

%-- STEP 3: Matix multiplication for LPU --%
var_lpu = J' * G * J;  % variance
s_lpu = sqrt(var_lpu);   % standard deviation
per_lpu = s_lpu / P_lpu * 100  % percent uncertainty


%== Monte Carlo ==========================================================%
%-- STEP 1: Generate samples --%
ns = 2e4;  % number of Monte Carlo samples
Nis = mvnrnd(Ni_bar', G, ns);  % input samples

Ps_mc = Nis(:,2) ./ Nis(:,1);  % STEP 2: Compute output for each sample

%-- STEP 3: Compute statistics --%
P_mc = mean(Ps_mc);
s_mc = std(Ps_mc);   % standard deviation
per_mc = s_mc / P_mc * 100  % percent uncertainty


%== Boostrapping =========================================================%
%-- STEP 1: Generate samples --%
ns = 2e4;  % number of Monte Carlo samples
r = randi(n, [size(Ni, 1), ns]);
idx = sub2ind(size(Ni), [ones(1,ns); 2 .* ones(1,ns)], r);
Nis = Ni(idx)';  % input samples

Ps_bs = Nis(:,2) ./ Nis(:,1);  % STEP 2: Compute output for each sample

%-- STEP 3: Compute statistics --%
P_bs = mean(Ps_bs);
s_bs = std(Ps_bs);   % standard deviation
per_bs = s_bs / P_bs * 100  % percent uncertainty


%== POST-PROCESSING ======================================================%
vec = linspace(-4 * s_lpu, 4 * s_lpu, 100) + P_lpu;
plot(vec, normpdf(vec, P_lpu, s_lpu) .* sqrt(2*pi) .* s_lpu, 'r');

vec2 = linspace(-6 * s_lpu, 6 * s_lpu, 75) + P_mc;
hold on;
[y, x] = histcounts(Ps_mc, vec2);
stairs([x(1), x], [0, y, 0] ./ max(y), 'k');
hold off;

vec2 = linspace(-6 * s_lpu, 6 * s_lpu, 75) + P_bs;
hold on;
[y, x] = histcounts(Ps_bs, vec2);
stairs([x(1), x], [0, y, 0] ./ max(y), 'g');
hold off;

xline(P_mc, 'k');
xline(P_lpu, 'r');
xline(P_bs, 'g');

xlim([-6 * s_lpu, 6 * s_lpu] + P_lpu);
xlabel('Penetration');
ylabel('Counts, probability density function');


