
% SDM2SMP  Calculate the GSD for the particle mass distribution from mobility distribution GSD.
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

function sm = sdm2smp(sd, prop)

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
sm = exp(log(sd) .* zet);

end

