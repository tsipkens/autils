
clear;
close all;
clc;

addpath cmap;


P = logspace(-3, 0, 200);

the = 1;
gam = 5;
Nvec = [1e3,1e4,1e5,1e6];


% Absolute uncertainties.
figure(1);
title('Absolute uncertainties');
xlabel('Penetration');
ylabel('std(P)');

cmap_sweep(length(Nvec), internet);
for ii=1:length(Nvec)
    loglog(P, P.* sqrt(the / Nvec(ii) .* (1 + 1 ./ P)));
    hold on;
end
hold off;

hold on;
for ii=1:length(Nvec)
    loglog(P, P.* sqrt(the / Nvec(ii) .* (1 + 1 ./ P) + ...
        (gam / Nvec(ii)) ^ 2 .* (1 + 1 ./ P .^ 2)));
end
hold off


% Relative uncertainties.
figure(2);
title('Relative uncertainties');
xlabel('Penetration');
ylabel('std(P)/P');

cmap_sweep(length(Nvec), internet);
for ii=1:length(Nvec)
    loglog(P, sqrt(the / Nvec(ii) .* (1 + 1 ./ P)));
    hold on;
end
hold off;

hold on;
for ii=1:length(Nvec)
    loglog(P, sqrt(the / Nvec(ii) .* (1 + 1 ./ P) + ...
        (gam / Nvec(ii)) ^ 2 .* (1 + 1 ./ P .^ 2)));
end
hold off
