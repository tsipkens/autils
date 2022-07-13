
% DVE2CHI  Compute the dynamic shape factor at a given mobility diameter.
%  
%  AUTHOR: Timothy Sipkens, 2022-07-13

function chi = dve2chi(dve, prop, f_iter)

%-- Parse inputs --------------------------%
if ~exist('f_iter', 'var'); f_iter = []; end
if isempty(f_iter); f_iter = 1; end  % default is iterative method

if ~isfield(prop, 'rhom')
    error('Shape factor calculation require rhom in prop.');
end
%-----------------------------------------%

dm = dve2dm(dve, prop);

chi = dm / dve;

if f_iter
    chi = (dm / dve * Cc(dve) / Cc(dm));
end

end
