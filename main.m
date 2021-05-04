

T = 298;
P = 101325;

prop.Dm = 3;  % mass-mobility exponent
prop.m0 = 4.7124e-25;  % mass-mobility pre-factor


[B, Zp] = dm2zp(100e-9, 1, T, P)

[B, Zp, d] = mp2zp(1e-18, 1, T, P, prop)

