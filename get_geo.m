

function [gmd, gsd, J] = get_geo(Ni, di, f_fit)

if ~exist('f_fit', 'var'); f_fit = []; end
if isempty(f_fit); f_fit = 0; end  % by default, use statistical definition

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

if f_fit  % instead, use lognormal fitting
    di = di(~isnan(Ni));  % truncate NaN values
    Ni = Ni(~isnan(Ni));

    x0 = [max(Ni) .* sqrt(2*pi) .* log(gsd), ...
        gmd, gsd];
    
    fun = @(x) Ni - ...
        normpdf(log(di), log(x(2)), log(x(3))) .* ...
        sqrt(2*pi) .* log(x(3)) .* x(1);
    
    opts.Display = 'none';
    [x1, ~, ~, ~, ~, ~, J] = lsqnonlin(fun, x0, [], [], opts);
    
    gmd = x1(2);
    gsd = x1(3);
end

end

