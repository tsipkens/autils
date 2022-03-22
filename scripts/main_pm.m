
clear;
close all;



dg = 75;
sg = 1.5;

chi = 1.0;


n = 220;
d = logspace(log10(10), log10(2e3), n)';

prop = massmob.gen('universal');


x0 = 1e6 .* normpdf(log(d), log(dg), log(sg)) .* (log(d(2)) - log(d(1)));

rng(0);
[xs, ~, Gx] = uq.add_noise(x0, 0.005, 1, 0.0001 .* max(x0), 4);
x = mean(x0, 2);
sx = sqrt(diag(Gx));
% sx = [sx; 0.0 .* prop.rho100; 0.0];

figure(1);
plot(d, xs, '.');
hold on;
plot(d, x, 'k');
hold off;
set(gca, 'XScale', 'log');


[M, Gm] = pm.pm_ni(x, d .* 1e-9, prop, sx);
pe = [sqrt(Gm) / M, M .* 1e12]  % 1e12 is mg

[M, Gm] = pm.pm_hc(x, d .* 1e-9, prop, sx);
pe = [sqrt(Gm) / M, M .* 1e12]



