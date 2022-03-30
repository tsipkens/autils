

clear;
close all;
clc;


P = 0.7;
sp = 3;


dg = 75;
sg = 1.65;

cmd = dg
mmd = hc(dg, sg, 3)



n = 2e3;
d = logspace(log10(10), log10(1000), n);

% Evaluate distribution.
x = normpdf(log(d), log(cmd), log(sg)) .* (sqrt(2*pi) * log(sg));
xm = normpdf(log(d), log(mmd), log(sg)) .* (sqrt(2*pi) * log(sg));


f = figure(1);
for ii=1:3
    subplot(1,3,ii);
    plot(d, x);
    hold on;
    plot(d, xm);
    hold off;
    
    ylim([0, 1.1]);
    xlim([min(d), max(d)]);
    set(gca, 'XScale', 'log');
    axis square;
end



% Pen A. 
pen = P .* ones(size(d));
subplot(1,3,1);
hold on;
plot(d, pen);
plot(d, pen .* x);
plot(d, pen .* xm);
hold off;



% Pen B.
dp = cmd .* 0.4;
pen = P .* normpdf(log(d), log(dp), log(sp)) .* (sqrt(2*pi) * log(sp));
subplot(1,3,2);
hold on;
plot(d, pen);
plot(d, pen .* x);
plot(d, pen .* xm);
hold off;



% Pen C.
dp = mmd .* 2.5;
pen = P .* normpdf(log(d), log(dp), log(sp)) .* (sqrt(2*pi) * log(sp));
subplot(1,3,3);
hold on;
plot(d, pen);
plot(d, pen .* x);
plot(d, pen .* xm);
hold off;


f.Position(3) = 1200;
