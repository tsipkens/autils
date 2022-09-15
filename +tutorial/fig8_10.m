
% A script demonstrating uncertainties in penetration.
%  Runtimes on the order of 1-2 minutes.
%  
%  AUTHOR: Timothy Sipkens: 2022-04-29

clear;
close all;
clc;
addpath cmap;

P = 0.9;  % penetration
N1 = 2e4;  % upstream number
N2 = P .* N1;  % downstream number

% Size of grid.
% n1 = 20;
% n2 = 21;
n1 = 105;
n2 = 120;

s1 = linspace(0.005, 0.4, n1);  % relative error, std(N1)/N1
s2 = linspace(0.005, 0.6, n2);  % relative error, std(N2)/N2



%%
% FIG 1. Linear error propagation (LUP) intervals.
sp_lup = sqrt(s1'.^2 + s2.^2);  % LUP
figure(1);
imagesc(s2, s1, sp_lup);
hold on;
contour(s2, s1, sp_lup, 0:0.1:1.1, 'w');
hold off;

% Format plot.
colorbar;
caxis([0, 1]);
addpath('cmap');
colormap(internet);
set(gca, 'YDir', 'normal');
ylabel('se(N_1)/N_1');
xlabel('se(N_2)/N_2');





%%
% Indices to save for map in FIG. 5. 
a1 = ceil(n1/2);  % 10
a2 = ceil(n2/2);  % 2

rng(1);  % set random number generator

% MAIN LOOP.
sp_prct = [];
sp_std = [];
ps = {};
smpl = {};
ns = 1e4;
for jj=1:n2
    for ii=1:n1
        [~, ~, smpl, ps] = uq.mc([N1, N2], ...
            [s1(ii) * N1, s2(jj) * N2], ...
            @(x) x(2,:) ./ x(1,:), 1e4, randi(1e5));
        
        ps2 = ps;
        ps2(any(smpl < eps)) = []; % remove negative number concentrations
        ps2 = rmoutliers(ps2, 'median', 'ThresholdFactor', 20);
        
        % Average percentile intervals to approx. std. dev.
        sp_prct(ii,jj) = (prctile(ps2, 83) - ...
            prctile(ps2, 17)) ./ P ./ 2;
        sp_std(ii,jj) = std(ps2);
        sk(ii,jj) = skewness(ps2);
        
        % Save samples for FIG. 5. 
        if and(ii == a2, jj == a1)
            smpla = smpl;
            psa = ps;
        end
    end
end

% FIG. 2. Monte Carlo intervals.
% Loop over std(N1)/N1 and std(N2)/N2.
figure(2);
imagesc(s2, s1, sp_std);
hold on;
contour(s2, s1, sp_std, 0:0.1:1.1, 'w');
hold off;

% Format plot.
colorbar;
caxis([0, 1]);
addpath('cmap');
colormap(internet);
set(gca, 'YDir', 'normal');
ylabel('se(N_1)/N_1');
xlabel('se(N_2)/N_2');


figure(3);
imagesc(s2, s1, sk);
hold on;
contour(s2, s1, sk, -0.4:0.4:5, 'k');
hold off;

% Format plot.
colorbar;
addpath('cmap');
colormap(flipud(internet));
set(gca, 'YDir', 'normal');
ylabel('se(N_1)/N_1');
xlabel('se(N_2)/N_2');


figure(4);
imagesc(s2, s1, (sp_std - sp_lup) ./ sp_lup);
hold on;
contour(s2, s1, (sp_std - sp_lup) ./ sp_lup, -1:0.1:1, 'k');
hold off;

% Format plot.
colorbar;
cmax = max(max((sp_std - sp_lup) ./ sp_lup));
clim([-cmax, cmax])
addpath('cmap');
colormap(curl);
set(gca, 'YDir', 'normal');
ylabel('se(N_1)/N_1');
xlabel('se(N_2)/N_2');


%%
%== FIG 5 ==%
% Plot mapping Nup to penetration.
l = length(smpla(1,:));
[y, idx] = sort(smpla(1,:) ./ mean(smpla(1,:)));

nc = 200;
cm = internet(nc);
ys = std(y);
ym = mean(y);
yc = (y - ym) ./ ys;
yc = min(3, max(-3, yc));
yc = ceil((0.5 + yc ./ 6) .* nc + eps);

figure(5);
subplot(1,5,2:4);
for ii=1:5:l
    plot([0;1], [y(ii); psa(1,idx(ii))], ...
        'Color', [cm(ceil(ii ./ l .* 200),:), 0.2], 'LineWidth', 1);
    hold on;
end
hold off;
xlim([0, 1]);
ylim([0, 2.2]);

% Plot left panel.
subplot(1,5,1);
[yh,xh] = histcounts(smpla(1,:) ./ mean(smpla(1,:)));
stairs([yh,0], xh, 'k');
set(gca, 'XDir', 'reverse','xtick',[]);
ylim([0, 2.2]);
xlabel('No. of samples');
ylabel('N_1 / mean(N_1)');

% Plot right panel.
subplot(1,5,5);
[yh,xh] = histcounts(psa(1,:));
stairs([yh,0], xh, 'k');
set(gca,'YAxisLocation','right','xtick',[]);
ylim([0, 2.2]);
xlabel('No. of samples');
ylabel('Penetration');


