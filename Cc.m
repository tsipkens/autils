
% CC  Function to evaluate Cunningham slip correction factor.
%  
%  CC = Cc(D) computes the Cunningham slip correction factor for
%  the provided mobility diameter, D. 
%  
%  CC = Cc(D,T,P) add inputs for the temperature in Kelvin, T, 
%  and pressure in atm., P. 
%  
%  ------------------------------------------------------------------------
% 
%  NOTE:
%   The temperature (T) and pressure (p) must both be specified 
%   before they are used. Otherwise, they are ignored and the 
%   code from Buckley et al. (2017) is used.
%  
%  AUTHOR: Timothy Sipkens, 2019-01-02

function Cc = Cc(d,T,p)

if nargin==1 % if P and T are not specified, use Buckley/Davies
    mfp = 66.5e-9; % mean free path
    
    % for air, from Davies (1945)
    A1 = 1.257;
    A2 = 0.4;
    A3 = 0.55;
    
else % from Olfert laboratory / Kim et al.
    S = 110.4; % temperature [K]
    mfp_0 = 6.730e-8; % mean free path of gas molecules in air [m]
    T_0 = 296.15; % reference temperature [K]
    p_0 = 101325; % reference pressure, [Pa] (760 mmHg to Pa)
    
    p = p*p_0;
    
    mfp = mfp_0*(T/T_0)^2*(p_0/p)*((T_0+S)/(T+S)); % mean free path
        % Kim et al. (2005) (doi:10.6028/jres.110.005), ISO 15900 Eqn 4
    
    A1 = 1.165;
    A2 = 0.483;
    A3 = 0.997/2;
    
end

Kn = (2*mfp)./d; % Knudsen number
Cc = 1 + Kn.*(A1 + A2.*exp(-(2*A3)./Kn)); % Cunningham slip correction factor

end
