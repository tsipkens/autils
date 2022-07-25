
% A script to generate Fig. 5 of Sipkens et al.

clear;
close all;
addpath cmap;

N = 1e4;
rho_100 = 800;
dg = 250;
sg = 1.6;

% Power law parameters.
prop = massmob.init('zet', 3, 'rho100', rho_100)

% Covariance information for power law.
var_rho = (0.005 .* prop.rho100) ^ 2 * 0;
var_zet = (0.005 .* prop.zet) ^ 2 * 0;
cov_rho_zet = 0. * sqrt(var_rho * var_zet);  % if correlated
G_rho_zet = [var_rho, cov_rho_zet; cov_rho_zet, var_zet];


n0 = 2e4;  % number of points in high diameter resolution distribution (for plotting)
d0 = logspace(-2, 6, n0)';  % vector of diameters (high res.)
dd0 = log(d0(2)) - log(d0(1));  % element width

n = 114;
d = logspace(log10(13.1), log10(763.5), n)';
dd = log(d(2)) - log(d(1));

p0 = N .* normpdf(log(d0), log(dg), log(sg)) .* dd0;
H0 = dm2mp(dg * 1e-9, prop) .* exp(prop.zet ^ 2 / 2 .* log(sg) .^ 2);
M0 = [sum(dm2mp(d0 .* 1e-9, prop) .* p0), N .* H0];
table(M0(1), M0(2), 'RowNames', {'M0'}, 'VariableNames', {'NI', 'SHC'})

r = randi(1e3);
p1 = N .* normpdf(log(d), log(dg), log(sg)) .* dd;
[p, Lp, outp] = uq.add_noise(p1, 0, 1, 1, 1, r);  % Poisson-Gaussian noise
Gp_inv = full(Lp' * Lp);
Gp = inv(Gp_inv);


di = d ./ 1e9;


%== NUMERICAL INTEGRATION (NI) ===========================================%
p_nim = N .* prop.m100 .* (d ./ 100) .^ prop.zet .* p;

[M_ni, s_ni, ~, J_ni] = pm.pm_ni(p, di, prop, ...
    blkdiag(Gp, G_rho_zet));
per_ni = s_ni / M_ni;

% (Ni, di, prop, G)
%== HATCH-CHOATE, SIMPLE (HCS) ===========================================%
[dg_hc, sg_hc] = get_geo(p, d);
p_hcn = N .* normpdf(log(d0), log(dg_hc), log(sg_hc)) .* dd;
p_hcm = N .* normpdf(log(d0), log(hc(dg_hc, sg_hc, prop.zet)), log(sg_hc)) .* dd;
[M_hcs, s_hcs, J_hcs] = pm.pm_hc(p, di, prop, ...
    blkdiag(Gp, G_rho_zet));
per_hcs = s_hcs / M_hcs;


figure(3);
dmax = 3e3;
stairs([d; d(end); dmax] - dd/2, [p ./ sum(p .* dd); 0; 0], 'k');
hold on;
stairs([d; 0; 0] - dd/2, [p_nim ./ sum(p .* dd) ./ M0(1); 0; 0], 'k');
plot(d0, p_hcn ./ sum(p_hcn .* dd0), 'Color', [0.1,0.2,0.8]);
plot(d0, p_hcm ./ sum(p_hcm .* dd0), 'Color', [0.1,0.2,0.8]);
plot(d0, p0 ./ sum(p0 .* dd0), 'r');
% stairs();
xline(max(d), 'k--');
hold off;

set(gca, 'XScale', 'log');
xlim([min(d), dmax]);


