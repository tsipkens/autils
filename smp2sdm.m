
% SMP2SDM  Calculate the GSD for the mobility distribution from particle mass distribution GSD.
%  
%  SM = sdm2smp(SD,PROP) converts the GSD using the mass-mobility
%  relation parameters in PROP (specifically PROP.zet). 
% 
%  SM = sdm2smp(SD,ZET) allows for explicit input of the mass-mobility
%  exponent, ZET. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function sd = smp2sdm(sm, prop)

% Parse inputs.
if isfloat(prop)
    zet = prop;
else
    if ~isfield(prop, 'zet')
        prop = massmob.init(prop); % else, make sure 'zet' is a field of prop.
    end
    zet = prop.zet;
end

% Use the mass-mobility relationship to get the new GSD.
sd = exp(log(sm) ./ zet);

end

