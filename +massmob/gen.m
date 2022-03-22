
% GEN  Fill in mass-mobility information using name-value pairs.
%  Includes computing prop.m0, which is used for mp2zp and mp2dm.
%  
%  NAME-VALUE options (mass-mobility exponent +1 other required):
%  
%   'm0' - mass of a 1 nm particle
%   'm100' - mass of a 100 nm particle
%   'rho0' - effective density of a 1 nm particle
%   'rho100' - effective density of a 100 nm particle
%   'zet' or 'Dm' - mass-mobility expeonent (required)
%  
%  NOTE: Can work in conjunction with other prop structures. Call 
%  `massmob.update(...)` to add mass-mobility parameters to an existing 
%  prop structure. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2021-03-25

function prop = gen(prop_str1, val1, str2, val2)

% Check for presents.
if ischar(prop_str1)
    switch prop_str1
        case {'NaCl', 'salt'}  % assume spheres with bulk density
            prop_str1 = 'zet';  val1 = 3;
            str2 = 'rho100';    val2 = 2160;
        case {'universal', 'soot'}  % universal soot relation (Olfert and Rogak)
            prop_str1 = 'zet';  val1 = 2.48;
            str2 = 'rho100';    val2 = 510;
    end
end

% Copy inputs to structure. 
if ischar(prop_str1)
    if ~exist('val2', 'var')
        error('Invalid input for mass-mobility relation.');
    else
        prop = struct();
        prop.(prop_str1) = val1;
        prop.(str2) = val2;
    end
else
    prop = prop_str1;
end

% Check afor mass-mobility exponent information.
if and(~isfield(prop, 'zet'), ~isfield(prop, 'Dm'))
    error('Mass-mobility exponent is required for mass-mobility relation.')
elseif ~isfield(prop, 'zet')
    prop.zet = prop.Dm;  % new standard in codes
else
    prop.Dm = prop.zet;  % for backwards compatibility
end

% Comptue properties. 
if ~isfield(prop, 'm0')
    if isfield(prop, 'm100')
        prop.m0 = prop.m100 * (1 / 100) ^ prop.zet;
        
    elseif isfield(prop, 'rho0')
        prop.m0 = prop.rho0 * pi / 6 * 1e-27;
        
    elseif isfield(prop, 'rho100')
        prop.m100 = prop.rho100 * pi / 6 * (100e-9) ^ 3;
        prop.m0 = prop.m100 * (1 / 100) ^ prop.zet;
        
    else
        error('Could not compute prop.m0.');
    end
end

prop.rho0 = prop.m0 * 6 / pi * 1e27;
prop.m100 = prop.m0 / (1 / 100) ^ prop.zet;
prop.rho100 = prop.m100 * 6 / pi / (100e-9 ^ 3);
prop.k = prop.m0;  % copy to k (alternative notation)

end
