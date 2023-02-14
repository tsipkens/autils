
% A script demonstrating error corruption of signals for lognormal distribution.
%  Roughly recreates Fig. 1.
%  
%  AUTHOR: Timothy Sipkens: 2022-06-23

clear;
close all;
addpath cmap;


sd = 90;  % seed for random number generation


%== DEFINE BASE SIGNAL/DISTRIBUTION ======================================%
% Distribution parameters.
N = 4e3;
dg = 75;
sg = 1.5;

% Discretization:
% Roughly similar to SMPS spacing.
n = 114;
d = logspace(log10(15), log10(400), n)';
dd = log(d(2)) - log(d(1));

p0 = N .* normpdf(log(d), log(dg), log(sg)) .* dd;
%=========================================================================%


%== POISSON NOISE SUBPLOT ================================================%
[p1, ~, G1] = uq.add_noise(p0, 0, 1, 0, 3, sd, 1);  % add Poisson noise

figure(1);
subplot(1,2,1);
plot(d, p0, 'k-');
cmap_sweep(3, flipud(internet));
hold on;
plot(d, p1, '.', 'MarkerSize', 8);
hold off;
xlim([min(d), max(d)]);
ylim([0, 200]);
set(gca, 'XScale', 'log');

%{
% Generate plot on right axis.
idx1 = 57;  idx2 = 36;   % indices in the signal
yvec = 0:0.01:200;  % for Gaussian representation
ypvec = 0:1:200;  % for Poisson representation

figure(2);
plot(normpdf(yvec, p0(idx1), sqrt(G1(idx1,idx1))), yvec, 'b-');
hold on;
plot(normpdf(yvec, p0(idx2), sqrt(G1(idx2,idx2))), yvec, 'r-');
ylim([0, 200]);

% Add Poisson representation.
ppdf = poisspdf(ypvec, p0(idx2));  % Compute Poisson distr.
stairs(ppdf, ypvec+0.5, 'r');
hold off;
%}


%== PGM ERRORS SUBPLOT ===================================================%
[p2, ~, G2] = uq.add_noise(p0, 0.25, 1, 0.2, 3, sd, 1);  % add PGM error
p2 = p2(:, [2,3,1]);  % re-order (just for ascending color in plot)

figure(1);
subplot(1,2,2);
plot(d, p0, 'k-');
cmap_sweep(3, flipud(internet));
hold on;
plot(d, p2, '.', 'MarkerSize', 8);
hold off;
xlim([min(d), max(d)]);
ylim([0, 200]);
set(gca, 'XScale', 'log');


