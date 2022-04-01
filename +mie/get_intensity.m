
% GET_INTENSITY  Gets the scattering intensity at a specific angle.
%  
%  AUTHOR: Timothy Sipkens, 29-03-2022

function [sn, sl, sr, the] = get_intensity(lam, d, n, nm, the)

%-- Parse inputs ------------------------------%
if ~exist('nm', 'var'); nm = []; end
if isempty(nm); nm = 1; end  % assume air as medium

if ~exist('the', 'var'); the = []; end
if isempty(the); the = linspace(0, pi, 300); end

lam = lam .* 1e-6;  % convert from m to micron
d = d .* 1e-6;
%----------------------------------------------%


%-- Pre-compute quantities --------------------%
m = n ./ nm;   % ratio of refractive indices
if all(size(m) == 1)
    if all(size(lam) == 1)  % if evaluating at a single wavelength
        m = m .* ones(size(d));
        lam = lam .* ones(size(d));
    else
        m = m .* ones(length(d), length(lam));
    end
end

x = pi .* d ./ (lam ./ nm);  % size parameter


%-- Main loop ---------------------------------%
s1 = zeros(size(x));  % initialize arrays
s2 = s1;

for ii=1:size(d, 1)  % loop over wavelength and diameter
    for jj=1:size(lam, 2)
        for kk=1:length(the)
            [s1(ii,jj,kk), s2(ii,jj,kk)] = mie.mie_s12(m(ii,jj), x(ii,jj), the(kk));
        end
    end
end

sl = squeeze(real(s1) .^ 2 + imag(s1) .^ 2);  % perpindicular
sr = squeeze(real(s2) .^ 2 + imag(s2) .^ 2);  % parallel
sn = (sr + sl) ./ 2;  % natural

end