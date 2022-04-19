

clear;
close all;
clc;



dg = 75;
sg = 1.86; %65;
% sg = 1.8;

chi = 1.0;

prop = massmob.init('salt');
rho = prop.rho100;

prop.chi = chi;



n = 2e3;
d = logspace(log10(10), log10(2e3), n);

f_aero = 0;

cmd = dg
[cmad, cmvd] = dm2da(cmd .* 1e-9, rho, chi, f_aero);
cmad = cmad .* 1e9
cmvd = cmvd .* 1e9

mmd = hc(dg, sg, 3)
smd = hc(dg, sg, 6)
[mmad, mmvd] = dm2da(mmd .* 1e-9, rho, chi, f_aero);

da = dm2da(d .* 1e-9, rho, chi, f_aero) .* 1e9;


sga = exp( ...
    (log(dm2da(exp(log(dg) .* 1.05) .* 1e-9, rho, chi, f_aero)) - ...
     log(dm2da(exp(log(dg) .* 0.95) .* 1e-9, rho, chi, f_aero))) ./ ...
    (0.1 .* log(dg)) .* ...
    log(sg))  % aerodynamic diameter, count distribution GSD

mmad = mmad .* 1e9
mmvd = mmvd .* 1e9


x = normpdf(log(d), log(cmd), log(sg)) .* (sqrt(2*pi) * log(sg));


xm = normpdf(log(d), log(mmd), log(sg)) .* (sqrt(2*pi) * log(sg));
xs = normpdf(log(d), log(smd), log(sg)) .* (sqrt(2*pi) * log(sg));

% Total scattering.
Qsca = mie.get_eff(532e-9, d' .* 1e-9, 1.5442 + 0j)';
Csca1 = Qsca .* (pi .* d .^ 2 ./ 4);

% 45 degree scattering.
Csca0 = mie.get_intensity(532e-9, d' .* 1e-9, 1.5442 + 0j, [], 45/180*pi)';

% 40-50 degree scattering.
Csca2 = mie.get_intensity(532e-9, d' .* 1e-9, 1.5442 + 0j, [], ...
    linspace(30, 60, 25)./180.*pi)';
Csca2 = mean(Csca2);

xs0 = x .* Csca0;
xs0 = xs0 ./ max(xs0);

xs1 = x .* Csca1;
xs1 = xs1 ./ max(xs1);

xs2 = x .* Csca2;
xs2 = xs2 ./ max(xs2);

figure(4);
plot(d, Csca0);
hold on;
plot(d, Csca1);
plot(d, d .^ 6 .* Csca0(1) ./ (d(1) .^ 6), 'k--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');


figure(1);
plot(d, x);
hold on;
plot(d, xm);
plot(d, xs);
plot(d, xs0);
plot(d, xs1);
plot(d, xs2);

xline(cmd);
xline(mmd);
xline(smd);

xline(exp(log(cmd) + log(sg)), 'g');
xline(exp(log(cmd) - log(sg)), 'g');

xline(cmad, 'r');
xline(mmad, 'r');
hold off;

set(gca, 'XScale', 'log');


%-{
% Figure used for translating log-scales.
figure(2);
d2 = [10:10:100, 200:100:1000, 2000:1000:4000];
da2 = dm2da(d2 .* 1e-9, rho, chi, f_aero) .* 1e9;
loglog(d2, da2, '.');
hold on;
loglog(mmd, mmad, 'r.');
loglog(cmd, cmad, 'g.');
hold off
xlim([min(d2), max(d2)]);
ylim([min(da2), max(da2)]);
%}


% Compare variants of aerodynamic diameter distributions.
% Second variant does not assume zet = 3.
figure(3);
plot(da, x);
hold on;
plot(da, normpdf(log(da), log(cmad), log(sga)) .* (sqrt(2*pi) * log(sga)));
hold off;
set(gca, 'XScale', 'log');


