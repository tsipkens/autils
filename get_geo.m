

function [gmd, gsd] = get_geo(Ni, di)

gmd = exp(sum(Ni ./ nansum(Ni) .* ...
        log(di), 'omitnan'));
gmd = real(gmd);

gsd = exp(sqrt(sum(Ni ./ nansum(Ni) .* ...
    (log(di) - log(gmd)) .^ 2, ...
    'omitnan')));
gsd = real(gsd);

% Cap lower limit of GSD, make monodisperse.
if gsd < exp(0.06)
    gsd = 1;
end

% if any(~isreal([gmd, gsd]))
%     disp('Imaginary GMD or GSD!');
% end

end

