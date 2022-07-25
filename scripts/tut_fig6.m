
% A script to generate Fig. 6 of Sipkens et al.

clear;
close all;
addpath cmap;

N = 2e4;
rho_100 = 2160;
dg = 75;
sg = 1.5; % 1.05

% Power law parameters.
prop = massmob.init('zet', 3, 'rho100', rho_100)

% Covariance information for power law.
var_rho = (0.005 .* prop.rho100) ^ 2;
var_zet = (0.005 .* prop.zet) ^ 2;
cov_rho_zet = 0 .* sqrt(var_rho * var_zet);  % if correlated
G_rho_zet = [var_rho, cov_rho_zet; cov_rho_zet, var_zet];


n0 = 2e4;  % overall number concentration
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
[M_ni, s_ni, ~, J_ni] = pm.pm_ni(p, di, prop, ...
    blkdiag(Gp, G_rho_zet));
per_ni = s_ni / M_ni;

%== HATCH-CHOATE, SIMPLE (HCS) ===========================================%
[M_hcs, s_hcs, ~, J_hcs] = pm.pm_hc(p, di, prop, ...
    blkdiag(Gp, G_rho_zet));
per_hcs = s_hcs / M_hcs;


%== HATCH-CHOATE, LOGNORMAL (HCL) ========================================%
[M_hcl, s_hcl] = pm.pm_hc_fit(p, di, prop, ...
    blkdiag(Gp, G_rho_zet));
per_hcl = s_hcl / M_hcl;


%== HATCH-CHOATE, SIMPLE, INDEPENDENT N (HCS2) ===========================%
[M_hcs_ind, s_hcs_ind] = pm.pm_hc_ind(p, di, prop, sum(p), ...
    blkdiag(Gp, G_rho_zet, sum(sum(Gp))));
per_hcs2 = s_hcs_ind / M_hcs_ind;


%== MONTE CARLO (NI) =====================================================%
x = [p1; prop.rho100; prop.zet];
Gx = blkdiag(Gp, G_rho_zet);

fk = @(x) x(end-1,:) .* (pi ./ 6 * (100e-9) ^ 3) .* (1 / 100) .^ x(end,:);
f = @(x) fk(x) .* sum((di .* 1e9) .^ x(end,:) .* x(1:(end-2), :));
[M_mcni, s_mcni, smpl, M_smpl] = uq.mc(x, Gx, f, 2e4, 6);
per_mcni = s_mcni / M_mcni;


%== MONTE CARLO (HCS) ====================================================
f_N = @(x) sum(x(1:(end-2),:));
f_dg = @(x) exp(1 ./ f_N(x) .* sum(x(1:(end-2),:) .* log(di)));
f_sg = @(x) sqrt(1 ./ f_N(x) .* ...
    sum(x(1:(end-2),:) .* (log(di) - log(f_dg(x))) .^ 2));
f_H = @(x) fk(x) .* (f_dg(x) .* 1e9) .^ x(end,:) .* exp(x(end,:) .^ 2 ./ 2 .* f_sg(x) .^ 2);
fh = @(x) f_N(x) .* f_H(x);
[M_mch, s_mch, smplh, M_smplh] = uq.mc(x, Gx, fh, 2e4, 6);
per_mch = s_mch / M_mch;



% Output results.
tbl = table([M_ni, M_hcs, M_hcs_ind, M_hcl, M_mcni, M_mch]' ./ M0(1), ...
    [s_ni, s_hcs, s_hcs_ind, s_hcl, s_mcni, s_mch]', ...
    100 .* [per_ni, per_hcs, per_hcs2, per_hcl, per_mcni, per_mch]', ...
    'RowNames', {'NI', 'HCS', 'HCS(IND)', 'HCL', 'MC(NI)', 'MC(HC)'}, ...
    'VariableNames', {'M/M0', 's', 'per'})


% Compare NI and MC(MI).
figure(4);
[hy, hx] = histcounts(M_smpl ./ M0(1));
[hyh] = histcounts(M_smplh ./ M0(1), hx);
vec = linspace(0, max(hx), 1e3);

stairs([hx(1), hx], [0, hy, 0], 'k');
hold on;
stairs([hx(1), hx], [0, hyh, 0], 'k');
plot(vec, normpdf(vec, M_mcni ./ M0(1), s_ni ./ M0(1)) .* sqrt(2*pi) .* s_ni ./ M0(1) .* max(hy));
plot(vec, normpdf(vec, M_mch ./ M0(1), s_hcs ./ M0(1)) .* sqrt(2*pi) .* s_ni ./ M0(1) .* max(hy));
plot(vec, normpdf(vec, 1, s_hcl ./ M0(1)) .* sqrt(2*pi) .* s_ni ./ M0(1) .* max(hy));
xline(1);
% xline((M_ni + s_ni) ./ M0(1), 'r');  % lines indicated standard devs.
% xline((M_hcs + s_hcs) ./ M0(1), 'g');
% xline((M_hcl + s_hcl) ./ M0(1), 'b');
% xline((M_ni - s_ni) ./ M0(1), 'r');
% xline((M_hcs - s_hcs) ./ M0(1), 'g');
% xline((M_hcl - s_hcl) ./ M0(1), 'b');
hold off;

xlim([0.7, 1.3]);


% Plot of scaled Jacobian for Ntot, NI, and HCS.
figure(1);
plot(d, ones(size(p1)) ./ N);
hold on;
plot(d, J_ni(1:end-2) ./ M0(1));
plot(d, J_hcs(1:end-2) ./ M0(2));
hold off;
set(gca, 'YScale', 'log', 'XScale', 'log');
xlim([min(d), max(d)]);
ylim([1e-5, 1e-2]);








