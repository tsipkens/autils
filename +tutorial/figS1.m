
% A script to generate Fig. 5 of Sipkens et al.

clear;
close all;
addpath cmap;

% Distribution parameters.
N = 1e4;  % Used to scale number distribution.
rho_100 = 800;  % Effective density at100 nm mobility diameter.
dg = 250;  % geometric mean diameter
sg = 1.6;  % goemetric standard deviation

% Mass-mobility parameters (display to console).
prop = massmob.init('zet', 3, 'rho100', rho_100)


n0 = 2e4;  % number of points in high diameter resolution distribution (for plotting)
d0 = logspace(-2, 6, n0)';  % vector of diameters (high res.)
dd0 = log(d0(2)) - log(d0(1));  % element width
Ni0 = N .* normpdf(log(d0), log(dg), log(sg)) .* dd0;  % "true" distribution for plotting


%== GENERATE SIGNALS =====================================================%
n = 114;
d = logspace(log10(13.1), log10(763.5), n)';
di = d ./ 1e9;
dd = log(d(2)) - log(d(1));

r = randi(1e3);
Ni1 = N .* normpdf(log(d), log(dg), log(sg)) .* dd;
[Ni, ~, Gp] = uq.add_noise(Ni1, 0, 1, 1, 1, r);  % Poisson-Gaussian noise
%=========================================================================%

% Compute M for noiseless data (to check that they are the same/scale below).
H0 = dm2mp(dg * 1e-9, prop) .* exp(prop.zet ^ 2 / 2 .* log(sg) .^ 2);
M0 = [sum(dm2mp(d0 .* 1e-9, prop) .* Ni0), N .* H0];  % [NI, SHC]
table(M0(1), M0(2), 'RowNames', {'M0'}, 'VariableNames', {'NI', 'SHC'})


%== STEP 0 ==%
% Covariance information for power law.
var_rho = (0.005 .* prop.rho100) ^ 2 * 0;
var_zet = (0.005 .* prop.zet) ^ 2 * 0;
cov_rho_zet = 0. * sqrt(var_rho * var_zet);  % if correlated
G_rho_zet = [var_rho, cov_rho_zet; cov_rho_zet, var_zet];

G = blkdiag(Gp, G_rho_zet);  % compile covariance for measurements and model params.


%== NUMERICAL INTEGRATION (NI) ===========================================%
%-- STEP 0: Expected values --%
[M_ni, s_ni, ~, J_ni] = pm.pm_ni(Ni, di, prop, G);  % using function instead
Ni_nim = N .* prop.m100 .* (d ./ 100) .^ prop.zet .* Ni;  % weighted distribution (for plotting)


%== HATCH-CHOATE, SIMPLE (HCS) ===========================================%
%-- STEP 0: Expected values --%
[dg_hc, sg_hc] = get_geo(Ni, d);
[M_hcs, s_hcs, ~, J_hcs] = pm.pm_hc(Ni, di, prop, G);
Ni_hcn = N .* normpdf(log(d0), log(dg_hc), log(sg_hc)) .* dd;  % weighted distribution (for plotting)
Ni_hcm = N .* normpdf(log(d0), log(hc(dg_hc, sg_hc, prop.zet)), log(sg_hc)) .* dd;


%== FIG ==%
figure(5);
dmax = 3e3;
stairs([d; d(end); dmax] - dd/2, [Ni ./ sum(Ni .* dd); 0; 0], 'k');
hold on;
stairs([d; 0; 0] - dd/2, [Ni_nim ./ sum(Ni .* dd) ./ M0(1); 0; 0], 'k');
plot(d0, Ni_hcn ./ sum(Ni_hcn .* dd0), 'Color', [0.1,0.2,0.8]);
plot(d0, Ni_hcm ./ sum(Ni_hcm .* dd0), 'Color', [0.1,0.2,0.8]);
plot(d0, Ni0 ./ sum(Ni0 .* dd0), 'r');
% stairs();
xline(max(d), 'k--');
hold off;

set(gca, 'XScale', 'log');
xlim([min(d), dmax]);


