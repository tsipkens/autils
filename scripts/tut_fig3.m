
% A script demonstrating error models.
%  
%  AUTHOR: Timothy Sipkens: 2022-04-19

clear;
close all;
clc;

% Load mu and sig from real SMPS signal.
load('data/smps1.mat');

[~, tau, the, gam] = uq.covp(mu1, 'pgm', sig1);


figure(2);

%-- Real SMPS signal represented with PGM error model --------------------%
plot(mu1, sig1 .^ 2, '.');
set(gca, 'XScale', 'log', 'YScale', 'log');

limy = ylim;
limx = xlim;
vec = logspace(log10(limx(1)), log10(limx(2)), 80);
hold on;
plot(vec, gam .^ 2 .* ones(size(vec)));  % Gaussian
plot(vec, the .* vec);  % Poisson
plot(vec, (tau .* vec) .^ 2);  % multiplicative
plot(vec, (tau .* vec) .^ 2 + the .* vec + gam .^ 2);  % PGM
hold off;

ylabel('Variance');
xlabel('Mean / Expected value')
ylim(limy);

%-- Simulation Poisson signal --------------------------------------------%
the2 = 10;
s2 = uq.add_noise(mu1 ./ 5, 0, the2, 0, 200);
hold on;
plot(mu1 ./ 5, var(s2, [], 2), '.');
plot(vec, the2 .* vec)
hold off;

