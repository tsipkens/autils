
% DVE2DM  Convert the mobility diameter to a volume-equivalent diameter. 
%  
%  DM = dve2dm(DVE,PROP) converts the mobility diameter, DM, to an 
%  volume equivalent diameter, DVE, using the material density and mass-
%  mobility information in PROP.
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2022-07-13

function dm = dve2dm(dve, prop)

dm = ((prop.rhom * pi / (6 * prop.k)) .* dve .^ 3) .^ (1 / prop.zet) .* 1e-9;

end
