
% USES FUNCTIONS FROM:
% Christian Maetzler, 
% MATLAB Functions for Mie Scattering and Absorption
% University of Bern, Institut fuer Angewandte Physik,
% Research Report No. 2002-08,June, 2002 
% http://arrc.ou.edu/~rockee/NRA_2007_website/Mie-scattering-Matlab.pdf

clear;
close all;
clc;

n = 1.5442 + 0j; % particle refractive index (salt)
% n = 2 + 1j; % particle refractive index (soot)

l = [0.532, 1.064] .* 1e-6; % [m] (starting wavelength, interval, ending wavelength)
d = logspace(log10(0.02), log10(2), 450)' .* 1e-6; % [m]  list them, [a b c], up to three particle diameters

[qsca, qabs] = mie.get_eff(l, d, n);


% Plot efficiencies.
figure(1);
clr = 'gr';

for jj=1:length(l)
    plot(d .* 1e6, qsca(:,jj), [clr(jj) '-'])
    hold on
    
    % Plot Rayleigh limit equivalent.
    plot(d .* 1e6, qsca(1,jj) ./ (d(1) .^ 4) .* d .^ 4, [clr(jj) '--'])
    
    %{
    plot(diadia, qabs(ii,:), [clr(ii) '-'])
    plot(diadia, qabs(ii,1) ./ diadia(1) .* diadia, [clr(ii) '--'])
    %}
end
hold off;
xlabel('Diameter d [\mum]')
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim([1e-4, 10]);


%%
% Angular scattering plots.
l2 = 523e-9;  % reduce to a single wavelength
d2 = logspace(log10(0.2), log10(2), 12)' .* 1e-6;

disp('Evaluate angular scattering...');
[in, ~, ~, the] = mie.get_intensity(l2, d2, n);
disp('DONE.');

login = log10(in);

figure(2);
cmap_sweep(length(d2), tokyo);
plot(the, in);
set(gca, 'YScale', 'log');
xlim([0, pi]);

