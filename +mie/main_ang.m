
% Computation and plot of Mie Power Scattering function for given 
% complex refractive-index ratio m=m'+im", size parameters x=k0*a, 
% according to Bohren and Huffman (1983) BEWI:TDD122
% C. Mätzler, May 2002.


clear;
close all;
clc;

n = 1.54;
nm = 1;
d = 0.4;
lam = 1;
x = pi .* d ./ (lam ./ nm);
nsteps = 300;

n=real(n); m2=imag(n);
nx=(1:nsteps); dteta=pi/(nsteps-1);
teta=(nx-1).*dteta;

for j = 1:nsteps
    
    u=cos(teta(j));
    [a(j),b(j)]=mie.mie_s12(n,x,acos(u));
    sl(j)= real(a(j)'*a(j));
    sr(j)= real(b(j)'*b(j));

end

yl=[teta, teta+pi; sl, fliplr(sl)]';  % perpindicular
yr=[teta, teta+pi; sr, fliplr(sr)]';  % parallel

% Natural.
sn = (sl + sr) ./ 2;
y2 = [teta teta+pi; sn, fliplr(sn)]';

[sn2, sl2, sr2, the2] = mie.get_intensity(lam, d, n, nm);

polar(yl(:,1),yl(:,2))
hold on;
polar(yr(:,1),yr(:,2))
polar(y2(:,1),y2(:,2))
polar(the2',squeeze(sn2), 'k--')
polar(the2',squeeze(sl2), 'k--')
polar(the2',squeeze(sr2), 'k--')
hold off;

title(sprintf('Mie angular scattering: m=%g+%gi, x=%g',n,m2,x));

xlabel('Scattering Angle')


