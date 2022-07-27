
% A script to demonstrate error propagation for total number concentration.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-09

clear;
close all;
addpath cmap;

%== SIGNAL GENERATION ====================================================%
% Can be replaced by your own signals.
% Parameters for distribution.
N0 = 1e5;  % total number (used to scale distribution and QoI)   < CHANGE THIS
dg = 100;  % GMD
sg = 1.6;  % GSD

% Generate diameter elements/discretization.
n = 114;  % number of diameters to have measurements at
d = logspace(log10(13.1), log10(763.5), n)';  % diameters
dd = log(d(2)) - log(d(1));  % element size in log-space

% Compute distribution and noisy data.
p = normpdf(log(d), log(dg), log(sg)) .* dd;  % noiseless particle size distribution
[Ni, ~, Gni0] = uq.add_noise(N0 .* p, 0, 1, 0, 20, 10);  % add noise
%=========================================================================%

%== STEP 0 ==%
% Get expected value/mean.
Ni_bar = mean(Ni, 2);

% Ways to compute covariance. 
Gni = Gni0;           % (1) uses true covariance from error model above
% Gni = cov(Ni');       % (2) use covariance computed from repeats
% Gni = uq.covp(Ni');   % (3) fit a PGM error model to repeats


% Compute the expected value for number concentration from noisy data.
N = sum(Ni_bar);  % compute number concentration
red = N ./ N0     % value relative to true quantity


%-- FIG: A single noisy signal used as input data --%
figure(1);
dmax = max(d);
stairs([d; d(end); dmax] - dd/2, [Ni(:,1) ./ sum(Ni(:,1) .* dd); 0; 0], 'k');
set(gca, 'XScale', 'log');
xlim([min(d), dmax]);


% Propagate errors via LPU.
J = ones(size(Ni_bar));  % STEP 2: Jacobian
Gn = J' * Gni * J;       % STEP 3: LPU matrix multiplication
rel = sqrt(Gn) ./ N0     % relative error


% Propagate errors via Monte Carlo.
Nis = mvnrnd(Ni_bar, Gni, 1e4)';  % sample number concentrations for each size bin
Ns = sum(Nis);  % samples of QoI


% Propagate errors via Monte Carlo, using a resampling method.
Nirs = uq.resample(Nis(:, 1:10), 1e4);  % resample using first few noisy signals
Nrs = sum(Nirs);  % samples of QoI


%-- FIG: Plot of uncertainties in N --%
figure(4);

% Plot Monte Carlo.
[hy, hx] = histcounts(Ns ./ N0);  % bin, relative to true N
dhx = hx(2) - hx(1);
ymax = max(hy);
hy = hy ./ max(hy);
stairs([hx, hx(end)] - dhx, [0, hy, 0], 'k');

% Plot bootstrapping.
[hy, hx] = histcounts(Nrs ./ N0);  % bin, relative to true N
dhx = hx(2) - hx(1);
hy = hy ./ ymax;
hold on;
stairs([hx, hx(end)] - dhx, [0, hy, 0], 'g');
hold off;

% Plot LPU-implied.
xvec = linspace(min(hx), max(hx), 150);
hold on;
plot(xvec, normpdf(xvec, N ./ N0, sqrt(Gn) ./ N0) .* sqrt(2*pi * Gn) ./ N0, 'r');
hold off;

ylabel('Counts, pdf');
xlabel('Total number concentration (scaled by true value), N/N0');
legend({'Monte Carlo (mean/var)', 'Bootstrapping', 'LPU-implied'});





%{
%== ADD LOOP OVER N0 ==%
figure(3);  % second figure
N0_vec = 10 .^ [4:1:6];
perj = [];  Gj = [];
for jj=1:length(N0_vec)
    N0 = N0_vec(jj);

    [Ni, ~, Gni] = uq.add_noise(N0 .* p, 0, 1, 0, 1, 10);  % add noise
    Nis = mvnrnd(Ni, Gni, 5e4)';  % sample number concentrations for each size bin
    Ns = sum(Nis);
    N = mean(Ns);
    
    Gn = J' * Gni * J;  % LPU

    % Propagate errors via Monte Carlo.
    Nis = mvnrnd(Ni, Gni, 2e5)';  % sample number concentrations for each size bin
    Ns = sum(Nis);
    
    if jj == 1
        [hy, hx] = histcounts(Ns ./ N0_vec(jj));  % bin, relative to true N
        dhx = hx(2) - hx(1);
        maxh = max(hy);
        maxGn = Gn;
        maxN = N0;
        xvec = linspace(min(hx), max(hx), 400);
    else
        hy = histcounts(Ns ./ N0, hx);
    end

    hy = hy ./ maxh;
    stairs([hx, hx(end)] - dhx, [0, hy, 0], 'k');
    
    hold on;
    plot(xvec, normpdf(xvec, N ./ N0, sqrt(Gn) ./ N0) .* sqrt(2*pi * maxGn) ./ maxN);
    xline((N + sqrt(Gn)) ./ N0, 'r');
    xline((N - sqrt(Gn)) ./ N0, 'r');
end
hold off;
%}
