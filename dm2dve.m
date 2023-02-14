
% DM2DVE  Convert the mobility diameter to a volume-equivalent diameter. 
%  
%  DM = dm2dve(DM,PROP) converts the mobility diameter, DM, to an 
%  volume equivalent diameter, DVE, using the material density and mass-
%  mobility information in PROP.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2022-07-13

function dve = dm2dve(dm, prop)

dve = dm .* (dm2rhoeff(dm, prop) ./ prop.rhom) .^ (1/3);

%{
% Alternate form.
b = (6/pi * prop.k * 1e9 ^ prop.zet / prop.rhom) ^ (1/3);
dve = b .* dm .^ (prop.zet/3);
%}

end
