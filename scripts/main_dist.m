

clear;
close all;
clc;



dg = 75;
sg = 1.65;

chi = 1.0;


n = 2e3;
d = logspace(log10(10), log10(1e3), n);



cmd = dg
[cmad, cmvd] = dm2da(cmd .* 1e-9, 2130, chi, 0);
cmad = cmad .* 1e9
cmvd = cmvd .* 1e9

mmd = hc(dg, sg, 3)
smd = hc(dg, sg, 6)
[mmad, mmvd] = dm2da(mmd .* 1e-9, 2130, chi, 0);

da = dm2da(d .* 1e-9, 2130, chi, 0) .* 1e9;



mmad = mmad .* 1e9
mmvd = mmvd .* 1e9


x = normpdf(log(d), log(cmd), log(sg)) .* (sqrt(2*pi) * log(sg));


xm = normpdf(log(d), log(mmd), log(sg)) .* (sqrt(2*pi) * log(sg));
xs = normpdf(log(d), log(smd), log(sg)) .* (sqrt(2*pi) * log(sg));


figure(1);
plot(d, x);
hold on;
plot(d, xm);
plot(d, xs);

xline(cmd);
xline(mmd);
xline(smd);

xline(cmad, 'r');
xline(mmad, 'r');
hold off;

set(gca, 'XScale', 'log');


figure(2);
d2 = [10:10:100, 200:100:1000];
da2 = dm2da(d2 .* 1e-9, 2130, chi, 0) .* 1e9;
loglog(d2, da2, '.');
xlim([min(d2), max(d2)]);
ylim([min(da2), max(da2)]);



