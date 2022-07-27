
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
% FIG. 2. Monte Carlo intervals.
% Loop over std(N1)/N1 and std(N2)/N2.
sp_prct = [];
sp_std = [];
ps = {};
smpl = {};
ns = 1e4;
for jj=1:n2
    for ii=1:n1
        [~, ~, smpl{ii,jj}, ps{ii,jj}] = uq.mc([N1, N2], ...
            [s1(ii) * N1, s2(jj) * N2], ...
            @(x) x(2,:) ./ x(1,:), 1e4);
        % [~, ~, smpl, ps] = uq.mc_pois([N1, N2], ...
        %     [s1(ii), s2(jj)], ...
        %     @(x) x(2,:) ./ x(1,:), 1e4);
        
        ps2 = ps;
        ps2{ii,jj}(any(smpl{ii,jj} < eps)) = []; % remove negative number concentrations
        ps2{ii,jj} = rmoutliers(ps2{ii,jj}, 'median', 'ThresholdFactor', 20);
        
        % Average percentile intervals to approx. std. dev.
        sp_prct(ii,jj) = (prctile(ps2{ii,jj}, 83) - prctile(ps2{ii,jj}, 17)) ./ P ./ 2;
        sp_std(ii,jj) = std(ps2{ii,jj});
        sk(ii,jj) = skewness(ps2{ii,jj});
    end
end

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


% Plot mapping Nup to penetration.
n1 = ceil(size(smpl,1)/2);
n2 = ceil(size(smpl,2)/2);
l = length(smpl{n1,n2}(1,:));
[y, idx] = sort(smpl{n1,n2}(1,:) ./ mean(smpl{n1,n2}(1,:)));

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
    plot([0;1], [y(ii); ps{n1,n2}(1,idx(ii))], ...
        'Color', [cm(ceil(ii ./ l .* 200),:), 0.2], 'LineWidth', 1);
    hold on;
end
hold off;
xlim([0, 1]);
ylim([0, 2.2]);

% Plot left panel.
subplot(1,5,1);
[yh,xh] = histcounts(smpl{end,end}(1,:) ./ mean(smpl{end,end}(1,:)));
stairs([yh,0], xh, 'k');
set(gca, 'XDir', 'reverse','xtick',[]);
ylim([0, 2.2]);
xlabel('No. of samples');
ylabel('N_1 / mean(N_1)');

% Plot right panel.
subplot(1,5,5);
[yh,xh] = histcounts(ps{n1,n2}(1,:));
stairs([yh,0], xh, 'k');
set(gca,'YAxisLocation','right','xtick',[]);
ylim([0, 2.2]);
xlabel('No. of samples');
ylabel('Penetration');


