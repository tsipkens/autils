

clear;
close all;
addpath cmap;



dg = 75;
sg = 1.8;

chi = 1.0;


n = 80;
d = logspace(log10(10), log10(2e3), n)';
d = logspace(log10(5), log10(2e3), 80)';

prop = massmob.init('salt');

% PFE parameters. 
npfe0 = 0.7
dmpps = 200;
spfe = 3.5;

x_up0 = 1e6 .* normpdf(log(d), log(dg), log(sg)) .* (log(d(2)) - log(d(1)));

Pi = normpdf(log(d), log(dmpps), log(spfe)) .* (sqrt(2*pi) * log(spfe));
Pi = Pi .* (1 - npfe0) ./ sum(Pi .* x_up0) .* sum(x_up0);  % normalize based on npfe0
x_down0 = x_up0 .* Pi;

mpfe0 = pfe.mpfe_hc(x_up0, x_down0, d, prop)  % noiseless MPFE


% Generate noisy signals.
% rng(0);
ns = 25;
ns2 = 3;
[xs, ~, Gx] = uq.add_noise([x_up0; x_down0], 0.005, 1, 0.0001 .* max(x_up0), ns, randi(1e6));
xs_up = xs(1:length(x_up0), :);
x_up = mean(xs_up(:, 1:ns2), 2); % x_up0; % mean(xs_up, 2);

xs_down = xs(length(x_down0)+1:end, :);
x_down = mean(xs_down(:, 1:ns2), 2);  % only consider first three in mean (more realistic)

% Standard deviations in counts.
sx_up = sqrt(diag(Gx(1:length(d), 1:length(d))));
sx_down = sqrt(diag(Gx(length(d)+1:end, length(d)+1:end)));


% Number-based PFE.
[npfe, s_npfe] = pfe.npfe(x_up, x_down, Gx ./ ns, []);
pe_npfe = [npfe, s_npfe / npfe, npfe - npfe0]

% MPFE via numerical integration.
[mpfe, s_mpfe] = pfe.mpfe_ni(x_up, x_down, d .* 1e-9, prop, Gx ./ ns, [], 0.01);
pe_mpfe_ni = [mpfe, s_mpfe / mpfe, mpfe - mpfe0]

% MPFE via Hatch-Choate.
[mpfe_hc, s_mpfe_hc] = pfe.mpfe_hc(x_up, x_down, d .* 1e-9, prop, Gx ./ ns, [], 0.01);
pe_mpfe_hc = [mpfe_hc, s_mpfe_hc / mpfe_hc, mpfe_hc - mpfe0]

% Scattering-based PFE via numerical integration.
ri = 1.5442 + 0j;  % particle refractive index (salt)
int = mie.get_intensity(532e-9, d .* 1e-9, ri, [], 45/180*pi);
scapfe_45 = pfe.scapfe_ni(x_up0, x_down0, int)


% Size-resolved filtration efficiencies.
% For averaged scans. 
[~, s_spfe0] = pfe.spfe(x_up0, x_down0, Gx ./ ns, []);  % noiseless uncertainties
[spfe, s_spfe] = pfe.spfe(x_up, x_down, Gx ./ ns, []);

% For individual scans.
[~, s_spfes0] = pfe.spfe(x_up0, x_down0, Gx, []);  % noiseless uncertainties
spfes = [];
for ii=1:size(xs_up, 2)
    spfes(:,ii) = pfe.spfe(xs_up(:,ii), xs_down(:,ii), Gx ./ ns, []);
end

% Plot size-resolved PFEs with uncertainty bounds.
figure(2);
cmap_sweep(size(xs_up, 2), internet);
for ii=1:size(spfes, 2)
    scatter(d, spfes(:, ii), 8, 'o', 'filled', 'MarkerFaceAlpha', 1);
    hold on;
end
scatter(d, spfe, 75, [0,0,0], 'o');
plot(d, 1 - Pi, 'k');
plot(d, 1 - Pi + 2 .* s_spfes0, 'k');
plot(d, 1 - Pi - 2 .* s_spfes0, 'k');
plot(d, 1 - Pi + 2 .* s_spfe0, 'k');
plot(d, 1 - Pi - 2 .* s_spfe0, 'k');
hold off;
set(gca, 'XScale', 'log');
ylim([0.5, 1.2]);
xlim([min(d), max(d)]);


figure(3);
plot(d, x_up0 ./ max(x_up0));
hold on;
plot(d, x_down0 ./ max(x_up0));
plot(d, Pi);
hold off;
set(gca, 'XScale', 'log');
xlim([min(d), max(d)]);


