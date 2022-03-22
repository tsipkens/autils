

clear;
close all;



dg = 75;
sg = 1.8;

chi = 1.0;


n = 220;
d = logspace(log10(10), log10(2e3), n)';

prop = massmob.gen('salt');

% PFE parameters. 
npfe0 = 0.8;
dmpps = 200;
spfe = 3;

x_up0 = 1e6 .* normpdf(log(d), log(dg), log(sg)) .* (log(d(2)) - log(d(1)));

Pi = normpdf(log(d), log(dmpps), log(spfe)) .* (sqrt(2*pi) * log(spfe));
Pi = Pi .* (1 - npfe0) ./ sum(Pi .* x_up0) .* sum(x_up0);  % normalize based on npfe0
x_down0 = x_up0 .* Pi;

mpfe0 = pfe.mpfe_hc(x_up0, x_down0, d, prop)  % noiseless MPFE


% Generate noisy signals.
% rng(0);
ns = 4;
[xs, ~, Gx] = uq.add_noise([x_up0; x_down0], 0.005, 1, 0.0001 .* max(x_up0), ns, randi(1e6));
xs_up = xs(1:length(x_up0), :);
x_up = mean(xs_up, 2); % x_up0; % mean(xs_up, 2);

xs_down = xs(length(x_down0)+1:end, :);
x_down = mean(xs_down, 2); % x_down0; % mean(xs_down, 2);

figure(1);
plot(d, xs_up, '.');
hold on;
plot(d, x_up, 'k');
plot(d, x_down, 'k');
plot(d, Pi .* max(x_up));
hold off;
set(gca, 'XScale', 'log');


[npfe, s_npfe] = pfe.npfe(x_up, x_down, Gx ./ ns, []);
pe = [s_npfe / npfe, npfe, npfe - npfe0]


[mpfe, s_mpfe] = pfe.mpfe_ni(x_up, x_down, d .* 1e-9, prop, Gx ./ ns, [], 0.);
pe = [s_mpfe / mpfe, mpfe, mpfe - mpfe0]

[mpfe_hc, s_mpfe_hc] = pfe.mpfe_hc(x_up, x_down, d .* 1e-9, prop, Gx ./ ns, [], 0.);
pe = [s_mpfe_hc / mpfe_hc, mpfe_hc, mpfe_hc - mpfe0]


% Size-resolved filtration efficiencies.
[spfe, s_spfe] = pfe.spfe(x_up, x_down, Gx ./ ns, []);

figure(2);
plot(d, 1 - Pi);
hold on;
plot(d, spfe, '.');
plot(d, 1 - Pi + 2 .* s_spfe, 'k');
plot(d, 1 - Pi - 2 .* s_spfe, 'k');
hold off;
set(gca, 'XScale', 'log');
ylim([0, 1.4]);


