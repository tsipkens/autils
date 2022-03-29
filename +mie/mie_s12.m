
% Computation of Mie Scattering functions S1 and S2
% for complex refractive index, m=m'+im", 
% size parameter, x=k0*a; and scattering angle, the;
% where k0=vacuum wave number, a=sphere radius;
% s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
% 
% C. Mätzler, May 2002

function [S1, S2] = mie_s12(m, x, the)

u = cos(the);

nmax = round(2+x+4*x^(1/3));

ab = mie.mie_abcd(m,x);
an = ab(1,:);
bn = ab(2,:);

pt = mie.mie_pt(u, nmax);
pin = pt(1,:);
tin = pt(2,:);

n = (1:nmax);
n2 = (2*n+1)./(n.*(n+1));
pin = n2.*pin;
tin = n2.*tin;
S1 = (an*pin'+bn*tin');
S2 = (an*tin'+bn*pin');

end