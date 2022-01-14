

T = 298;
P = 1;

% prop = prop_massmob('Dm', 3, 'rho100', 900);
prop = prop_massmob('Dm', 2.48, 'rho100', 510);  % universal relation


[B, Zp, d] = mp2zp(1e-18, 1, T, P, prop)

[B, Zp] = dm2zp(d, 1, T, P)

[da, dve] = dm2da(d, 2130, 1.12, [], T, P)  % for salt
