
% A script to test most of the functions.
% Intended to check for errors following refactoring, not accuracy.

clear;
clc;

prop = massmob.init('universal')

massmob.add(prop, 'zet', prop.zet, 'rho100', prop.rho100 * 2)
massmob.update(prop, 'zet', 2.4)

dg = 100;
sg = 1.5;

da = dm2da(dg * 1e-9, 2160) * 1e9
dm = da2dm(da * 1e-9, 2160) * 1e9

mp = dm2mp(dg * 1e-9, prop) * 1e18
dm = mp2dm(mp * 1e-18, prop) * 1e9

[B, Zp] = dm2zp(dm, 1)

rho = rhoeff(dm * 1e-9, mp * 1e-18)

dq = hc(dg, sg, 3)
N = hc(dg, sg, 0, 1)  % integrates to unity, as expected
d1 = hc(dg, sg, 1, 1)  % mean, the same as `exp(log(dg) + log(1.5) ^ 2 / 2)`
Mq = hc(dg, sg, prop.zet, prop.k)


% Test functions requiring distributions.
d = logspace(log10(dg) - 3 .* log10(sg), ...
             log10(dg) + 3 .* log10(sg), 600)';  % diameters
dd = log(d(2)) - log(d(1));
p = normpdf(log(d), log(dg), log(sg)) .* N .* dd;

Mni = pm.pm_ni(p, d .* 1e-9, prop)
Mhc = pm.pm_hc(p, d .* 1e-9, prop)
Mhcl = pm.pm_hc_fit(p, d .* 1e-9, prop)

pfe0 = 0.8;
spfe = pfe.spfe(p, (1 - pfe0) .* p);
npfe = pfe.npfe(p, (1 - pfe0) .* p)
mpfe_ni = pfe.mpfe_ni(p, (1 - pfe0) .* p, d, prop)
mpfe_hc = pfe.mpfe_hc(p, (1 - pfe0) .* p, d, prop)
scapfe_hc = pfe.scapfe_ni(p, (1 - pfe0) .* p, ones(size(p)))



