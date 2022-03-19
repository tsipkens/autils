
% A script to test most of the functions.

clear;
clc;

prop = massmob.gen('universal')

massmob.update(prop, 'zet', prop.zet, 'rho100', prop.rho100 * 2)

dg = 100;

da = dm2da(dg * 1e-9, 2160) * 1e9
dm = da2dm(da * 1e-9, 2160) * 1e9

mp = dm2mp(dg * 1e-9, prop) * 1e18
dm = mp2dm(mp * 1e-18, prop) * 1e9

[B, Zp] = dm2zp(dm, 1)

rho = rho_eff(dm * 1e-9, mp * 1e-18)

dq = hc(dg, 1.5, 3)
N = hc(dg, 1.5, 0, 1)  % integrates to unity, as expected
d1 = hc(dg, 1.5, 1, 1)  % mean, the same as `exp(log(dg) + log(1.5) ^ 2 / 2)`
Mq = hc(dg, 1.5, 3, prop.zet)
