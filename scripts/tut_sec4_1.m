
% A script to demonstrate error propagation for total number concentration.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-09

clear;
close all;
addpath cmap;

% Parameters for distribution.
N0 = 1e4;  % total number (used to scale distribution and QoI)
dg = 100;  % GMD
sg = 1.6;  % GSD

% Generate diameter elements/discretization.
n = 114;  % number of diameters to have measurements at
d = logspace(log10(13.1), log10(763.5), n)';  % diameters
dd = log(d(2)) - log(d(1));  % element size in log-space

% Compute distribution and noisy data.
p = normpdf(log(d), log(dg), log(sg)) .* dd;  % noiseless particle size distribution
[Ni, ~, Gni] = uq.add_noise(N0 .* p, 0, 1, 0, 1, 10);  % add noise

% Compute number concentration from noisy data.
N = sum(Ni);  % compute number concentration
red = N ./ N0  % value relative to true quantity

%-- FIG: Noisy input data --%
figure(1);
dmax = max(d);
stairs([d; d(end); dmax] - dd/2, [Ni ./ sum(Ni .* dd); 0; 0], 'k');

set(gca, 'XScale', 'log');
xlim([min(d), dmax]);



% Propagate errors via LPU.
J = ones(size(Ni));  % Jacobian
Gn = J' * Gni * J;  % LPU
rel = sqrt(Gn) ./ N0  % relative error

% Propagate errors via Monte Carlo.
Nis = mvnrnd(Ni, Gni, 1e4)';  % sample number concentrations for each size bin
Ns = sum(Nis);

%-- FIG: Plot of uncertainties in N --%
figure(3);
[hy, hx] = histcounts(Ns ./ N0);  % bin, relative to true N
dhx = hx(2) - hx(1);
hy = hy ./ max(hy);
stairs([hx, hx(end)] - dhx, [0, hy, 0], 'k');

xvec = linspace(min(hx), max(hx), 150);
hold on;
plot(xvec, normpdf(xvec, N ./ N0, sqrt(Gn) ./ N0) .* sqrt(2*pi * Gn) ./ N0);
hold off;





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
