
% UPDATE  Update mass-mobility parameters in prop structure.
%  Uses a single name-value pair, holding another name-value pair constant.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-24

function prop = update(prop, f, v, fc)

%-- Parse inputs ------------------------%
if ~exist('fc', 'var'); fc = []; end
if isempty(fc); fc = 'rho100'; end  % use rho100 as constant field

%-- Update structure --------------------%
if or(strcmp(f, 'Dm'), strcmp(f, 'zet'))  % update exponent
    prop = massmob.init(fc, prop.(fc), 'zet', v);
    
else  % update pre-factor
    prop = massmob.init(f, v, 'zet', prop.zet);
    
end

end
