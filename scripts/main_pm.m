
clear;
close all;



dg = 75;
sg = 1.5;

chi = 1.0;


n = 220;
d = logspace(log10(10), log10(2e3), n)';

prop = massmob.init('universal');


x0 = 1e6 .* normpdf(log(d), log(dg), log(sg)) .* (log(d(2)) - log(d(1)));

% rng(0);
ns = 4;
[xs, ~, Gx] = uq.add_noise(x0, 0.005, 1, 0.0001 .* max(x0), ns, randi(1e6));
x = mean(xs, 2);
sx = sqrt(diag(Gx) ./ ns);
% sx = [sx; 0.0 .* prop.rho100; 0.0];

M0 = sum(x) * hc(dg / 100, sg, prop.zet, pi / 6 * prop.rho100 * (100^3) * 1e-27) * 1e12

figure(1);
plot(d, xs, '.');
hold on;
plot(d, x, 'k');
hold off;
set(gca, 'XScale', 'log');


[M, sm] = pm.pm_ni(x, d .* 1e-9, prop, sx);
pe = [sm / M, M .* 1e12, M .* 1e12 - M0]  % 1e12 is mg

[M, sm] = pm.pm_hc(x, d .* 1e-9, prop, sx);
pe = [sm / M, M .* 1e12, M .* 1e12 - M0]



