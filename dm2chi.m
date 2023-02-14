
% DM2CHI  Compute the dynamic shape factor at a given mobility diameter.
%  
%  AUTHOR: Timothy Sipkens, 2022-07-13

function chi = dm2chi(dm, prop, f_iter, varargin)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end  % default is iterative method

if ~isfield(prop, 'rhom')
    error('Shape factor calculation requires rhom in prop.');
end
%-----------------------------------------%

dve = dm2dve(dm, prop);

chi = dm ./ dve;
if f_iter
    chi = chi ./ Cc(dm, varargin{:}) .* Cc(dve, varargin{:});
end

%{
% Alternate form.
b = (6/pi * prop.k * 1e9 ^ prop.zet / prop.rhom) ^ (1/3);
chi = dm .^ (1 - prop.zet/3) ./ b;
if f_iter
    chi = chi ./ Cc(dm) .* Cc(b .* dm .^ (prop.zet/3));
end
%}

end
