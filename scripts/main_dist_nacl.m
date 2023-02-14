

clear;
close all;
clc;



dg = 75;
sg = 1.86; %65;
% sg = 1.8;

prop = massmob.init('salt');
prop.rhom = 2160;
% Differs from web app based on mass-mobility parameters.


n = 2e3;
d = logspace(log10(10), log10(2e3), n);

f_iter = 1;

cmd = dg
cmad = dm2da(cmd .* 1e-9, prop, f_iter);
cmvd = dm2dve(cmd .* 1e-9, prop);
cmad = cmad .* 1e9
cmvd = cmvd .* 1e9

mmd = hc(dg, sg, 3)
smd = hc(dg, sg, 6)
mmad = dm2da(mmd .* 1e-9, prop, f_iter);
mmvd = dm2dve(mmd .* 1e-9, prop);

da = dm2da(d .* 1e-9, prop, f_iter) .* 1e9;


sga = sdm2sda(sg, dg * 1e-9, prop, f_iter)  % aerodynamic diameter, count distribution GSD

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
da2 = dm2da(d2 .* 1e-9, prop, f_iter) .* 1e9;
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


