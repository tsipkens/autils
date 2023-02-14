
% CC  Function to evaluate Cunningham slip correction factor.
%  
%  CC = Cc(D) computes the Cunningham slip correction factor for
%  the provided mobility diameter, D, in nm.  By default, uses 'davies'.
%  
%  CC = Cc(D,T,P) add inputs for the temperature, T, in Kelvin and
%  pressure, P, in atm. By default, uses 'kim'.
%  
%  CC = Cc(D,OPT) uses the string in OPT to specify the slip correction
%  parameters.
%  
%  CC = Cc(D,T,P,OPT) uses the string in OPT to specify the slip 
%  correction parameters in addition to using T and P.
%  
%  ------------------------------------------------------------------------
% 
%  NOTE:
%   The temperature (T) and pressure (p) must both be specified 
%   before they are used. Otherwise, they are ignored and the 
%   code from Buckley et al. (2017) is used.
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function Cc = Cc(d, T_opt, p, opt)

if ~exist('T_opt', 'var'); T_opt = []; end

if ischar(T_opt); opt = T_opt;
elseif ~isempty(T_opt); T = T_opt; end

if ~exist('opt','var'); opt = []; end
if isempty(opt)
    if ~exist('p', 'var')
        opt = 'Davies';
    else
        opt = 'Kim';
    end
end
opt = lower(opt);


switch opt

    case 'davies'
    % For air from Davies (1945) as adopted by Buckley et al. (2017).
    % Default if P and T are not specified.

        mfp = 66.5e-9;  % mean free path [nm]
        A1 = 1.257;
        A2 = 0.4;
        A3 = 0.55;


    case {'hinds', 'allen', 'raabe', 'allen-raabe'}
    % For air from Allen and Raabe (1982, 1985) and cited in Hinds.
    % Paramters are adjusted given that computation below is in terms of
    % Knuden number (i.e., A1 is 2.34/2 and A2 is 1.05/2). 

        mfp = 65.1e-9;  % mean free path [nm]
        A1 = 1.17;
        A2 = 0.525;
        A3 = 0.39;
        

    case {'kim', 'iso'}
    % For air from Kim et al./ISO 15900 as adapted from Olfert et al.
    % Default if P and T are specified.

        S = 110.4;       % Sutherland constant [K]
        mfp_0 = 67.3e-9; % mean free path of gas molecules in air [m]
        T_0 = 296.15;    % reference temperature [K]
        p_0 = 101325;    % reference pressure, [Pa] (760 mmHg to Pa)
        
        p = p * p_0;
        
        % Kim et al. (2005) (doi:10.6028/jres.110.005), ISO 15900 Eqn 4
        % Correct default mean free path.
        mfp = mfp_0 * (T/T_0)^2 * (p_0/p) * ((T_0 + S)/(T + S));
        
        A1 = 1.165;
        A2 = 0.483;
        A3 = 0.997/2;
end

Kn = (2 * mfp) ./ d;  % Knudsen number
Cc = 1 + Kn .* (A1 + A2 .* exp(-(2 * A3) ./ Kn));  % Cunningham slip correction factor

end
