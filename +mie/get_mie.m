
% GET_MIE  Evaluates the sca/abs/ext efficiencies using Mie theory.
%  
%  QSCA = mie.get_mie(L, D, N) computes the scattering efficiency as a
%  function fo wavelengths, L; diameters, D; and particle refractive 
%  indices, N. The refractice indices can have imaginary components.
%  
%  [QSCA, QABS, QEXT] = mie.get_mie(...) adds outputs for the absorption
%  and extinction efficiencies.
%  
%  [...] = mie.get_mie(L, D, N, NM, FV)
%  fv            = volume fraction of spheres in medium (eg., fv = 0.05)
%  lambda        = wavelength in um (eg., lambda = 0.633)
%  dia           = sphere diameter in um (eg., dia_um = 0.0500)
%  np            = particle refractive index (eg. polystyrene = 1.57)
%  nm            = medium refractive index (eg., water = 1.33)
%                  Note: npar and nmed can be imaginary numbers.
%  
%  USES:
%   mie.m, which uses mie_abcd.m, from Maetzler (2002)
%  
%  AUTHORS:
%   Timothy Sipkens, 2022 (modified)
%   Steven Jacques, 2009 (original)

function [qsca, qabs, qext] = get_mie(lam, d, n, nm, fv)

%-- Parse inputs ------------------------------%
if ~exist('nm', 'var'); nm = []; end
if isempty(nm); nm = 1; end  % assume air as medium

if ~exist('fv', 'var'); fv = []; end
if ~isempty(fv)
    Vsphere = (4/3*pi) .* (d./2) .^ 3;  % volume of sphere
    rho = fv ./ Vsphere;  % #/um^3, concentration of spheres (only used in debugging)
    disp('Max. rho = ');
    disp(max(rho));
end

lam = lam .* 1e-6;  % convert from m to micron
d = d .* 1e-6;
%----------------------------------------------%


%-- Pre-compute quantities --------------------%

m = n ./ nm;   % ratio of refractive indices
if all(size(m) == 1)
    m = m .* ones(length(d), length(lam));
end

x = pi .* d ./ (lam ./ nm);  % size parameter


%-- Main loop ---------------------------------%
qsca = zeros(size(x));  % initialize arrays
qabs = qsca;  qext = qsca;  g = qsca;

for ii=1:length(d)  % loop over wavelength and diameter
    for jj=1:length(lam)
        out = mie.mie(m(ii,jj), x(ii,jj))';  % <----- Matlzer's subroutine
        % u = [real(m) imag(m) x qext qsca qabs qb asy qratio];
        
        %-- Unpack output -------%
        qext(ii,jj) = out(4);  % extinction efficiency, Qext
        qsca(ii,jj) = out(5);  % scattering efficiency, Qsca
        g(ii,jj)    = out(8);  % anisotropy, g
        qabs(ii,jj) = out(6);  % absorption efficiency, Qabs
    end
end

%{
% Other parameters, not currently used.
A       = pi*d^2/4;           % geometrical cross-sectional area, um^2
sigma_s = qsca*A;               % scattering cross-section, um^2
mus     = sigma_s*rho*1e4;      % scattering coeff. cm^-1
musp    = mus*(1-g);            % reduced scattering coeff. cm^-1

musgp = real([mus g musp]);
%}

end

