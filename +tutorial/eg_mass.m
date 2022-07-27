
% Example of computing uncertainties for total mass concentration.

clear;
close all;

% Mass-mobility relation parameters.
rho100 = 800;  % effective density at 100 nm mobility diameter.
zet = 3;  % mass-mobiltiy exponent
k = pi / 6 * rho100 * (100^3) * 1e-27;


% Load Poisson-Gaussian corrupted signals. 
% Contains Ni, d, dd, dg, sg, and N.
load('+tutorial/data/sim_pois_gaus');
H0 = k .* (dg / 100) .^ zet .* ...
    exp(1/2 .* zet .^ 2 .* (log(sg) .^ 2));  % "true" Hatch-Choate factor
M0 = N .* H0;  % true mass (as simulated signals, not possible for experiments)


%== STEP 0: Covariances and expected values ==============================%
Ni_bar = mean(Ni, 2);  % get mean of measurements
G_Ni = cov(Ni');  % get covariance

% Covariance information for power law.
var_rho = (0.005 .* rho100) ^ 2 * 0;
var_zet = (0.005 .* zet) ^ 2 * 0;
cov_rho_zet = 0. * sqrt(var_rho * var_zet);  % if correlated
G_rho_zet = [var_rho, cov_rho_zet; cov_rho_zet, var_zet];

G = blkdiag(G_Ni, G_rho_zet);  % compile covariance for measurements and model params.


%== NUMERICAL INTEGRATION (NI) ===========================================%
%-- STEP 0: Expected values --%
M_ni = k .* sum((d ./ 100) .^ zet .* Ni_bar);  % expected mass

%-- STEP 2: Jacobian --%
J_ni = [k .* (d ./ 100) .^ zet; ...  % for measurements
    M_ni ./ rho100; ...          % for mass-mobility exponent, zet
    k .* sum((d ./ 100) .^ zet .* Ni_bar .* log(d ./ 100))];  % for pre-factor, k

%-- STEP 3: Matix multiplication for LPU --%
var_ni = J_ni' * G * J_ni;  % variance
s_ni = sqrt(var_ni);   % standard deviation


%== HATCH-CHOATE, SIMPLE (HCS) ===========================================%
%-- STEP 0: Expected values --%
[dg_hc, sg_hc] = get_geo(Ni_bar, d);
H = k .* (dg_hc / 100) .^ zet .* ...
    exp(1/2 .* zet .^ 2 .* (log(sg_hc) .^ 2));  % expected Hatch-Choate factor
M_hcs = H .* sum(Ni_bar);  % expected mass

%-- STEP 2: Jacobian --%
Di = 1 + zet .* (log(d) - log(dg_hc)) + ...
    zet ^ 2 ./ 2 .* ((log(d) - log(dg_hc)) .^ 2 + log(sg_hc) .^ 2);
J_hcs = [H .* Di; ...
    M_hcs ./ rho100;  ...
    M_hcs * (zet * log(sg_hc) ^ 2 + log(dg ./ 100))];

%-- STEP 3: Matix multiplication for LPU --%
var_hcs = J_hcs' * G * J_hcs;  % variance
s_hcs = sqrt(var_hcs);   % standard deviation


%== POST-PROCESSING ======================================================%
% Table showing results.
per_ni = s_ni / M_ni;  % percent uncertainties
per_hcs = s_hcs / M_hcs;
table([M0; nan; nan], [M_ni; s_ni; per_ni * 100], [M_hcs; s_hcs; per_hcs * 100], ...
    'RowNames', {'M', 'Std', 'Std [%]'}, 'VariableNames', {'True', 'NI', 'SHC'})

figure(1);
plot(d, Ni, '.');
hold on;
plot(d, Ni_bar, 'k', 'LineWidth', 1);
hold off;
set(gca, 'XScale', 'log');
xlabel('Mobility diameter [nm]');
ylabel('Number concentration [counts/cc]');

figure(2);
s_max4 = 4 * max(s_ni, s_hcs);
vec = linspace(-s_max4, s_max4, 100) + M0;  % for x-axis of plots
plot(vec, normpdf(vec, M_ni, s_ni), 'r');
hold on;
plot(vec, normpdf(vec, M_hcs, s_hcs), 'b');
hold off;

xline(M0, 'k--');
xline(M_ni, 'r');
xline(M_hcs,'b');

xlabel('Total particulate mass, M');
ylabel('Probability density function (pdf)');
legend({'Numerical integration', 'Hatch-Choate', ...
    'M_0', 'M_{ni}', 'M_{hcs}'});


