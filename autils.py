
import numpy as np

from scipy.optimize import fsolve, minimize, root, least_squares
from scipy.stats import norm

from autils import props


def textdone():
    print('\r' +'\033[32m' + '^ DONE!' + '\033[0m' + '\n')


#== Functions for mass-mobility relations ======================#
def massmob_init(prop_str1=None, val1=None, str2=None, val2=None, d=100, rhom=None):
    """
    Fill in mass-mobility information using name-value pairs.
    Includes computing prop['m0'], which is used for mp2zp and mp2dm.

    NAME-VALUE options (mass-mobility exponent + 1 other required):
    - 'm0': mass of a 1 nm particle
    - 'm100': mass of a 100 nm particle
    - 'rho0': effective density of a 1 nm particle
    - 'rho100': effective density of a 100 nm particle
    - 'zet' or 'Dm': mass-mobility exponent (required)

    PROP = init(..., d) adds an input for diameter to be used with
    STR = 'rhod' or STR = 'md', corresponding to effective density or mass at a diameter of D.

    PROP = init(..., d, rhom) adds an input for the material density.
    
    AUTHOR: Timothy Sipkens, 2021-03-25
    """
    
    # Check for presets
    if isinstance(prop_str1, str):
        if prop_str1 in ['NaCl', 'salt']:  # assume spheres with bulk density
            prop_str1 = 'zet'
            val1 = 3
            str2 = 'rho100'
            val2 = 2160
        elif prop_str1 in ['universal', 'soot']:  # universal soot relation (Olfert and Rogak)
            prop_str1 = 'zet'
            val1 = 2.48
            str2 = 'rho100'
            val2 = 510
            rhom = 1860  # added at the bottom
        elif prop_str1 == 'water':  # water spheres
            prop_str1 = 'zet'
            val1 = 3
            str2 = 'rho100'
            val2 = 1000
            rhom = 1000  # added at the bottom
        elif prop_str1 == 'santovac':  # water spheres
            prop_str1 = 'zet'
            val1 = 3
            str2 = 'rho100'
            val2 = 1198
            rhom = 1198  # added at the bottom

    # Copy inputs to structure
    if isinstance(prop_str1, str):
        if val2 is None:
            raise ValueError('Invalid input for mass-mobility relation.')
        else:
            prop = {}
            prop[prop_str1] = val1
            prop[str2] = val2
    else:
        prop = prop_str1

    # Check for mass-mobility exponent information
    if 'zet' not in prop and 'Dm' not in prop:
        raise ValueError('Mass-mobility exponent is required for mass-mobility relation.')
    elif 'zet' not in prop:
        prop['zet'] = prop['Dm']  # new standard in codes
    else:
        prop['Dm'] = prop['zet']  # for backwards compatibility

    # Build structure
    if 'm0' not in prop:
        if 'm100' in prop:
            prop['m0'] = prop['m100'] * (1 / 100) ** prop['zet']
        elif 'md' in prop:
            prop['m0'] = prop['md'] * (1 / (d * 1e9)) ** prop['zet']
        elif 'rho0' in prop:
            prop['m0'] = prop['rho0'] * 3.14159 / 6 * 1e-27
        elif 'rho100' in prop:
            prop['m100'] = prop['rho100'] * 3.14159 / 6 * (100e-9) ** 3
            prop['m0'] = prop['m100'] * (1 / 100) ** prop['zet']
        elif 'rhod' in prop:
            prop['md'] = prop['rhod'] * 3.14159 / 6 * d ** 3
            prop['m0'] = prop['md'] * (1 / (d * 1e9)) ** prop['zet']
        else:
            raise ValueError('Could not compute prop["m0"].')

    # Fill out parameters
    prop['m100'] = prop['m0'] / (1 / 100) ** prop['zet']
    prop['rho0'] = prop['m0'] * 6 / 3.14159 * 1e27
    prop['rho100'] = prop['m100'] * 6 / 3.14159 / (100e-9 ** 3)
    prop['k'] = prop['m0']  # copy to k (alternative notation)

    # If an additional diameter was specified, reproduce
    if d != 100:
        prop['d'] = d
        prop['md'] = prop['m0'] / (1 / (d * 1e9)) ** prop['zet']
        prop['rhod'] = prop['md'] * 6 / 3.14159 / (d ** 3)

    # Check for rhom (material density) from parsing preset and add
    if rhom is not None:
        prop['rhom'] = rhom

    return prop


def massmob_add(prop, f1, v1=None, f2=None, v2=None):
    """
    Add or update mass-mobility parameters in an existing prop dictionary.

    This function can add mass-mobility information to a dictionary containing 
    flow rates, temperature, pressure, and other parameters for evaluating 
    the transfer functions of classifiers.

    Parameters:
    prop (dict): The prop structure (dictionary) to update.
    f1 (str): Name of particle type or field to add/update.
    v1 (any): Value corresponding to the field f1. If f1 is a string of particle type, this can be None.
    f2 (str): Second field to add/update (if using name-value pairs).
    v2 (any): Value corresponding to the field f2.

    Returns:
    dict: Updated prop dictionary with mass-mobility parameters.
    """

    # Remove the relevant fields to be replaced
    fields = ['zet', 'Dm', 'm0', 'k', 'rho0', 'm100', 'rho100']
    for field in fields:
        if field in prop:
            del prop[field]

    if v1 is None:  # If second argument is string of particle type
        prop1 = massmob_init(f1)  # Assuming init() is defined elsewhere
        prop['zet'] = prop1['zet']
        prop['rho100'] = prop1['rho100']

    else:  # If name-value pairs are provided
        # Add new values
        prop[f1] = v1
        prop[f2] = v2

    # Fill out remaining values
    prop = massmob_init(prop)

    return prop


def massmob_update(prop, f, v, fc=None):
    """
    Update mass-mobility parameters in the prop structure.
    Uses a single name-value pair, holding another name-value pair constant.
    
    Args:
    - prop (dict): The properties dictionary to update.
    - f (str): The field to update.
    - v (float): The new value for the field.
    - fc (str, optional): The field to keep constant. Defaults to 'rho100'.
    
    Returns:
    - dict: Updated properties dictionary.
    """
    if fc is None:
        fc = 'rho100'  # default to 'rho100' if not provided

    if f in ['Dm', 'zet']:  # update exponent
        prop = massmob_init(fc, prop[fc], 'zet', v)
    else:  # update pre-factor
        prop = massmob_init(f, v, 'zet', prop['zet'])
    
    return prop


#== Slip correction ============================================#
def cc(d, T_opt=None, p=None, opt=None):
    """
    Compute the Cunningham slip correction factor for the provided mobility diameter, d, in nm.
    
    Parameters:
        d (float or numpy array): Mobility diameter in nm.
        T_opt (float, optional): Temperature in Kelvin. If not provided, uses default settings.
        p (float, optional): Pressure in atm. If not provided, uses default settings.
        opt (str, optional): String to specify the slip correction parameters ('davies', 'kim', etc.).
    
    Returns:
        Cc (float or numpy array): Cunningham slip correction factor.
    """

    # Handle optional arguments and defaults
    if T_opt is None:
        T_opt = None
    
    if isinstance(T_opt, str):
        opt = T_opt
    elif T_opt is not None:
        T = T_opt

    if opt is None:
        if p is None:
            opt = 'davies'
        else:
            opt = 'kim'
    
    opt = opt.lower()

    # Set parameters based on the selected model
    if opt == 'davies':
        # For air from Davies (1945) as adopted by Buckley et al. (2017)
        mfp = 66.5e-9  # Mean free path [nm]
        A1 = 1.257
        A2 = 0.4
        A3 = 0.55

    elif opt in ['hinds', 'allen', 'raabe', 'allen-raabe']:
        # For air from Allen and Raabe (1982, 1985) and cited in Hinds
        mfp = 65.1e-9  # Mean free path [nm]
        A1 = 1.17
        A2 = 0.525
        A3 = 0.39

    elif opt in ['kim', 'iso']:
        # For air from Kim et al./ISO 15900 as adapted from Olfert et al.
        S = 110.4       # Sutherland constant [K]
        mfp_0 = 67.3e-9 # Mean free path of gas molecules in air [m]
        T_0 = 296.15    # Reference temperature [K]
        p_0 = 101325    # Reference pressure, [Pa]

        p = p * p_0  # Convert pressure from atm to Pa
        
        # Kim et al. (2005) (doi:10.6028/jres.110.005), ISO 15900 Eqn 4
        # Correct default mean free path.
        mfp = mfp_0 * (T / T_0)**2 * (p_0 / p) * ((T_0 + S) / (T + S))
        
        A1 = 1.165
        A2 = 0.483
        A3 = 0.997 / 2

    # Compute the Knudsen number
    Kn = (2 * mfp) / d
    
    # Compute the Cunningham slip correction factor
    Cc = 1 + Kn * (A1 + A2 * np.exp(-(2 * A3) / Kn))

    return Cc

#== Aerodynamic diameter conversions ===========================#
def da2dm(da, prop, f_iter=1, *args):
    """
    Convert aerodynamic diameter to mobility diameter.
    
    Parameters:
        da (float or numpy array): Aerodynamic diameter in meters.
        prop (dict): A dictionary containing the particle properties (e.g., density).
        f_iter (int, optional): Flag indicating whether to use simple (0) or iterative (1) evaluation. Default is 1 (iterative).
        *args: Additional arguments passed to the slip correction function Cc if needed.

    Returns:
        dm (float or numpy array): Mobility diameter in meters.
    """

    # Density of water in kg/m^3
    rho0 = 1e3

    # Function to calculate the effective density using dm2rhoeff (assumed to be implemented)
    def fun_simple(dm):
        return 1e9 * (dm * 1e-9 * np.sqrt(dm2rhoeff(dm * 1e-9, prop) / rho0) - da)

    # Solve for mobility diameter without iteration (simple method)
    dm = fsolve(fun_simple, da * 1e9) * 1e-9

    # Iterative method (more precise)
    if f_iter:
        def fun_iter(dm):
            return 1e9 * (dm * 1e-9 * np.sqrt(dm2rhoeff(dm * 1e-9, prop) / rho0 *
                           cc(dm * 1e-9, *args) / cc(da, *args)) - da)
        
        dm = fsolve(fun_iter, dm * 1e9) * 1e-9

    # Ensure result is real
    dm = np.real(dm)

    return dm

def da_rhoeff2dm(da, rho, f_iter=1, *args):
    """
    Convert aerodynamic diameter and effective density to mobility diameter.
    
    Parameters:
        da (float or numpy array): Aerodynamic diameter in meters.
        rho (float or numpy array): Effective density of same size as da or scalar.
        f_iter (int, optional): Flag indicating whether to use simple (0) or iterative (1) evaluation. Default is 1 (iterative).
        *args: Additional arguments passed to the slip correction function Cc if needed.

    Returns:
        dm (float or numpy array): Mobility diameter in meters.
    """

    # Density of water in kg/m^3
    rho0 = 1e3

    # Calculate the mobility diameter directly without iteration (simple method)
    dm = da * np.sqrt(rho0 / rho)

    # Iterative method (more precise)
    if f_iter:
        def fun_iter(dm):
            return 1e9 * (dm * 1e-9 * np.sqrt(rho / rho0 *
                           cc(dm * 1e-9, *args) / cc(da, *args)) - da)
        
        dm = fsolve(fun_iter, dm * 1e9) * 1e-9

    # Ensure result is real
    dm = np.real(dm)

    return dm

def da2dve(da, prop, f_iter=1):
    """
    Convert the aerodynamic diameter to a volume-equivalent diameter.
    
    Parameters:
    da (float or ndarray): Aerodynamic diameter.
    prop (dict): Dictionary containing the material density and other properties.
    f_iter (int, optional): Flag indicating whether to use the iterative method (default: 1).
    
    Returns:
    dve (float or ndarray): Volume-equivalent diameter.
    """
    rho0 = 1e3  # Density of water
    
    if f_iter:
        # Iterative method
        def fun(dve):
            return (dve * np.sqrt(prop['rhom'] / rho0 / dve2chi(dve, prop, f_iter) * cc(dve) / cc(da)) - da) * 1e9
    else:
        # Simple method
        def fun(dve):
            return (da / np.sqrt(prop['rhom'] / rho0 / dve2chi(dve, prop, f_iter)) - dve) * 1e9
    
    # Solve for dve using fsolve
    dve = fsolve(fun, da)
    
    return dve


def rhoeff(d, m):
    """
    Computes the effective density from mass and mobility diameter.
    
    Parameters:
        d: float or array-like
            Mobility diameter in meters.
        m: float or array-like
            Particle mass in kilograms.
    
    Returns:
        rho: float or array-like
            Effective density in kg/m^3.
    """
    return 6 * m / (np.pi * d**3)


#== Mobility diameter conversions ==============================#
def dm2chi(dm, prop, f_iter=1, *args):
    """
    Compute the dynamic shape factor at a given mobility diameter.

    Parameters:
    dm: Mobility diameter
    prop: Dictionary with properties including 'rhom'
    f_iter: Flag to determine if iterative method is used (default is iterative method)
    *args: Additional arguments for Cc function (if required)

    Returns:
    chi: Dynamic shape factor
    """
    # Check if 'rhom' is in the property dictionary
    if 'rhom' not in prop:
        raise ValueError('Shape factor calculation requires "rhom" in prop.')

    # Compute volume-equivalent diameter
    dve = dm2dve(dm, prop)

    # Initial chi calculation
    chi = dm / dve
    if f_iter:
        chi = chi / cc(dm, *args) * cc(dve, *args)

    # Alternative form (commented in MATLAB code)
    # b = (6 / np.pi * prop['k'] * 1e9 ** prop['zet'] / prop['rhom']) ** (1 / 3)
    # chi = dm ** (1 - prop['zet'] / 3) / b
    # if f_iter:
    #     chi = chi / cc(dm) * cc(b * dm ** (prop['zet'] / 3))

    return chi


def dm2da(dm, prop, f_iter=1, *args):
    """
    Convert the mobility diameter (dm) to an aerodynamic diameter (da).

    Parameters:
    dm (float or ndarray): Mobility diameter.
    prop (dict): Dictionary containing properties for conversion.
    f_iter (int, optional): Flag indicating whether to use simple (0) or iterative (1) evaluation. Default is 1.
    args: Additional arguments for the Cunningham slip correction factor (Cc).

    Returns:
    da (float or ndarray): Aerodynamic diameter.
    """
    rho0 = 1e3  # density of water in kg/m^3

    # Compute effective density
    rhoeff = dm2rhoeff(dm, prop)  # effective density from dm
    da = dm * np.sqrt(rhoeff / rho0)  # simple aerodynamic diameter

    # Alternate, iterative method
    if f_iter:
        def func(da1):
            return 1e9 * (da * np.sqrt(cc(dm, *args) / cc(da1, *args)) - da1)

        # Solve for aerodynamic diameter using fsolve
        da = fsolve(func, da)

    return da


def dm2dve(dm, prop):
    """
    Convert the mobility diameter to a volume-equivalent diameter.
    
    Parameters:
    dm (float or ndarray): Mobility diameter in meters.
    prop (dict): Dictionary containing the material density and mass-mobility properties.
    
    Returns:
    dve (float or ndarray): Volume-equivalent diameter in meters.
    """
    # Compute the volume-equivalent diameter
    return dm * (dm2rhoeff(dm, prop) / prop['rhom']) ** (1/3)


def dm2mp(d, prop):
    """
    Computes the effective density from mass and mobility diameter.
    
    Parameters:
        d: float or array-like
            Mobility diameter in meters.
        prop: dict
            Properties containing relevant physical parameters.
    
    Returns:
        m: float or array-like
            Particle mass in kilograms.
    """
    if not (type(prop) == dict):
        prop = props.massmob_init(prop)

    return prop['m0'] * (d * 1e9) ** prop['zet']


def dm2rhoeff(d, prop):
    """
    Get effective density from mobility diameter using mass-mobility relation.
    
    Parameters:
        d: float or array-like
            Mobility diameter (in meters).
        prop: dict
            Properties containing relevant physical parameters.
    
    Returns:
        rho: float or array-like
            Effective density.
    """
    return rhoeff(d, dm2mp(d, prop))


def dm2zp(dm, z=1, T=None, p=None):
    """
    Calculate electric mobility from mobility diameter and charge state.

    Parameters:
    dm (float or ndarray): Mobility diameter in meters.
    z (int, optional): Integer charge state. Defaults to 1.
    T (float, optional): Temperature in Kelvin. Defaults to None.
    p (float, optional): Pressure in atm. Defaults to None.

    Returns:
    B (float or ndarray): Mechanical mobility.
    Zp (float or ndarray): Electromobility.
    """
    # Define constants
    e = 1.6022e-19  # electron charge [C]

    # Default gas viscosity for Buckley/Davies
    if T is None or p is None:
        mu = 1.82e-5  # gas viscosity [Pa*s]
        B = cc(dm) / (3 * np.pi * mu * dm)  # mechanical mobility
    else:
        # Constants for Olfert laboratory / Kim et al.
        S = 110.4  # temperature [K]
        T_0 = 296.15  # reference temperature [K]
        vis_23 = 1.83245e-5  # reference viscosity [kg/(m*s)]
        mu = vis_23 * ((T / T_0) ** 1.5) * ((T_0 + S) / (T + S))  # gas viscosity
        B = cc(dm, T, p) / (3 * np.pi * mu * dm)  # mechanical mobility

    # Calculate electromobility
    Zp = B * e * z

    return B, Zp


def dm_da2mp(dm, da, *varargs):
    """
    Compute mass from mobility and aerodynamic diameters.

    Parameters:
    dm (float or ndarray): Mobility diameter in meters.
    da (float or ndarray): Aerodynamic diameter in meters.
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    mp (float or ndarray): Particle mass in kg.
    """
    return da**2 * dm * (np.pi * 1e3 / 6 * cc(da, *varargs) / cc(dm, *varargs))


def dm_da2rhoeff(dm, da, *varargs):
    """
    Compute effective density from mobility and aerodynamic diameters.

    Parameters:
    dm (float or ndarray): Mobility diameter in meters.
    da (float or ndarray): Aerodynamic diameter in meters.
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    rho (float or ndarray): Effective density in kg/m^3.
    """
    # Compute effective density
    return 1000 * (da / dm)**2 * cc(da, *varargs) / cc(dm, *varargs)


def dm_mp2da(dm, mp, f_iter=1, *varargs):
    """
    Compute aerodynamic diameter from mobility diameter and mass.

    Parameters:
    dm (float or ndarray): Mobility diameter in meters.
    mp (float or ndarray): Particle mass in kilograms.
    f_iter (int, optional): Flag indicating whether to use iterative method.
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    da (float or ndarray): Aerodynamic diameter in meters.
    """

    # Compute simple volume-equivalent and aerodynamic diameters.
    da = dm**(-1/2) * mp**(1/2) * np.sqrt(6 / (np.pi * 1e3))

    if f_iter:
        # Define the function to solve
        def fun(da1):
            return 1e9 * (da * np.sqrt(cc(dm, *varargs) / cc(da1, *varargs)) - da1)

        # Solve for aerodynamic diameter
        da = fsolve(fun, da)

    return da


#== Volume-equivalent diameter conversions =======================#
def dve2chi(dve, prop, f_iter=1):
    """
    Compute the dynamic shape factor at a given volume-equivalent diameter (dve).
    
    Parameters:
    dve (float or ndarray): Volume-equivalent diameter.
    prop (dict): Dictionary containing properties like 'rhom'.
    f_iter (int, optional): Flag indicating whether to use the iterative method (default: 1).
    
    Returns:
    chi (float or ndarray): Dynamic shape factor.
    """
    # Check if 'rhom' exists in the prop dictionary
    if 'rhom' not in prop:
        raise ValueError('Shape factor calculation requires rhom in prop.')

    # Compute the mobility diameter from the volume-equivalent diameter
    dm = dve2dm(dve, prop)

    # Compute the dynamic shape factor
    chi = dm / dve

    # If iterative method is requested, refine the calculation
    if f_iter:
        chi = (dm / dve) * (cc(dve) / cc(dm))

    return chi


def dve2da(dve, prop, f_iter=1, *varargs):
    """
    Convert volume-equivalent diameter to aerodynamic diameter.

    Parameters:
    dve (float or ndarray): Volume-equivalent diameter in meters.
    prop (dict): Dictionary containing properties including 'rhom'.
    f_iter (int, optional): Flag indicating whether to use iterative method.
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    da (float or ndarray): Aerodynamic diameter in meters.
    """
    rho0 = 1e3  # Density of water in kg/m^3

    # Compute simple volume-equivalent and aerodynamic diameters
    chi = dve2chi(dve, prop, f_iter)
    da = dve * np.sqrt(prop['rhom'] / rho0 / chi)

    if f_iter:
        # Define the function to solve
        def fun(da1):
            return 1e9 * (dve * np.sqrt(prop['rhom'] / rho0 / chi * 
                                         cc(dve, *varargs) / cc(da1, *varargs)) - da1)

        # Solve for aerodynamic diameter
        da = fsolve(fun, da)

    return da


def dve2dm(dve, prop):
    """
    Convert the volume-equivalent diameter (dve) to a mobility diameter (dm).

    Parameters:
    dve (float or ndarray): Volume-equivalent diameter.
    prop (dict): Dictionary containing properties like 'rhom', 'k', and 'zet'.
    
    Returns:
    dm (float or ndarray): Mobility diameter.
    """
    # Check if 'rhom' exists in the prop dictionary
    if 'rhom' not in prop:
        raise ValueError('The dve2dm function requires prop["rhom"].')

    # Compute the mobility diameter from the volume-equivalent diameter
    dm = ((prop['rhom'] * np.pi / (6 * prop['k'])) * dve ** 3) ** (1 / prop['zet']) * 1e-9

    return dm


#== Single particle mass conversions =============================#
def mp2dm(m, prop):
    """
    Calculate mobility diameter from particle mass using mass-mobility relation.

    Parameters:
    m (float or ndarray): Particle mass in kilograms.
    prop (dict): Dictionary containing mass-mobility properties with 'm0' and 'zet'.

    Returns:
    float or ndarray: Mobility diameter in meters.
    """
    if 'm0' not in prop:
        prop = massmob_init(prop)

    # Compute mobility diameter and return
    return 1e-9 * (m / prop['m0']) ** (1 / prop['zet'])


def mp_da2dm(mp, da, f_iter=1, *varargs):
    """
    Compute mobility diameter from particle mass and aerodynamic diameter.

    Parameters:
    mp (float or ndarray): Particle mass in kilograms.
    da (float or ndarray): Aerodynamic diameter in meters.
    f_iter (int, optional): Flag indicating whether to use iterative method.
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    dm (float or ndarray): Mobility diameter in meters.
    """
    # Compute simple volume-equivalent and aerodynamic diameters
    dm = da**(-2) * mp * (6 / (np.pi * 1e3))

    if f_iter:
        # Define the function to solve
        def fun(dm1):
            return 1e9 * (dm * (cc(dm1, *varargs) / cc(da, *varargs)) - dm1)

        # Solve for mobility diameter
        dm = fsolve(fun, dm)

    return dm


def mp2zp(m, z, T=None, P=None, prop=None):
    """
    Calculate electric mobility from particle mass.

    Parameters:
    m (float or ndarray): Particle mass in kilograms.
    z (int or ndarray): Integer charge state.
    T (float, optional): Temperature in Kelvin.
    P (float, optional): Pressure in atm.
    prop (dict): Dictionary containing mass-mobility properties with 'm0' and 'zet'.

    Returns:
    tuple: Mechanical mobility (B), electromobility (Zp), and mobility diameter (d).
    """
    if prop is None or 'm0' not in prop:
        raise ValueError("prop must contain 'm0'.")

    # Compute mobility diameter
    d = mp2dm(m, prop)

    # Compute mechanical and electrical mobility
    if T is None or P is None:
        B, Zp = dm2zp(d, z)
    else:
        B, Zp = dm2zp(d, z, T, P)

    return B, Zp, d


#== Mobility > mobility diameter ================================#
def zp2dm(Zp, z, T=None, P=None):
    """
    Calculate mobility diameter from a vector of mobilities.
    
    Parameters:
    Zp : array-like
        The electric mobility values.
    z : int
        Integer charge state.
    T : float, optional
        Temperature in Kelvin. If not specified, defaults to using Buckley/Davies method.
    P : float, optional
        Pressure in atm. If not specified, defaults to using Buckley/Davies method.
        
    Returns:
    dm : array-like
        The mobility diameters.
    B : array-like
        Mechanical mobility.
    """
    e = 1.6022e-19  # Electron charge [C]
    B = Zp / (e * z)  # Mechanical mobility
    
    if T is None or P is None:
        # If P and T are not specified, use Buckley/Davies method
        mu = 1.82e-5  # Gas viscosity [Pa*s]
        Cc_fun = lambda dm: cc(dm)
    else:
        # If P and T are specified
        S = 110.4  # Temperature [K]
        T_0 = 296.15  # Reference temperature [K]
        vis_23 = 1.83245e-5  # Reference viscosity [kg/(m*s)]
        mu = vis_23 * ((T / T_0) ** 1.5) * ((T_0 + S) / (T + S))  # Gas viscosity

        Cc_fun = lambda dm: cc(dm, T, P)

    
    dm0 = 1 / (3 * np.pi * mu * B)  # Initial estimate
    dm0 = dm0 * Cc_fun(dm0)

    def objective(dm):
        return 1e9 * np.linalg.norm(B - Cc_fun(np.abs(dm)) / (3 * np.pi * mu * np.abs(dm)))
    
    result = minimize(objective, dm0, method='Nelder-Mead', options={'disp': False})
    
    return result.x, B


#== Distribution manipulation ===================================#
def sdm2sda(sd, d, prop, f_iter=0):
    """
    Calculate the GSD for the aerodynamic distribution from mobility distribution GSD.

    Parameters:
    sd (float or ndarray): GSD for the mobility distribution.
    d (float or ndarray): Mobility diameter.
    prop (dict): Dictionary containing mass-mobility properties.
    f_iter (int, optional): Flag for iterative correction (default is 0).

    Returns:
    float or ndarray: GSD for the aerodynamic distribution.
    """
    d1 = np.exp(np.log(d) * 1.01)
    d2 = np.exp(np.log(d) * 0.99)

    da1 = dm2da(d1, prop, f_iter)
    da2 = dm2da(d2, prop, f_iter)

    sa = np.exp(
        (np.log(da1) - np.log(da2)) /
        (0.02 * np.log(d)) * np.log(sd)
    )

    return sa


def sdm2smp(sd, prop):
    """
    Calculate the GSD for the particle mass distribution from mobility distribution GSD.

    Parameters:
    sd (float or ndarray): GSD for the mobility distribution.
    prop (dict or float): Dictionary containing mass-mobility properties or the mass-mobility exponent directly.

    Returns:
    float or ndarray: GSD for the particle mass distribution.
    """
    if isinstance(prop, (int, float)):  # Check if prop is a float or int (i.e., the exponent zeta)
        zet = prop
    else:
        if 'zet' not in prop:
            prop = massmob_init(prop)  # Ensure 'zet' is a field of prop
        zet = prop['zet']

    # Calculate the new GSD using the mass-mobility relationship
    return np.exp(np.log(sd) * zet)


def smp2sdm(sm, prop):
    """
    Calculate the GSD for the mobility distribution from particle mass distribution GSD.
    
    Parameters:
    sm : array-like
        The GSD for the particle mass distribution.
    prop : dict or float
        If a dictionary, should contain the key 'zet'. If a float, is the mass-mobility exponent ZET.
        
    Returns:
    sd : array-like
        The GSD for the mobility distribution.
    """
    if isinstance(prop, float):
        zet = prop
    else:
        if 'zet' not in prop:
            prop = massmob_init(prop)  # Ensure 'zet' is a field of prop
        zet = prop['zet']
    
    # Use the mass-mobility relationship to get the new GSD.
    return np.exp(np.log(sm) / zet)
    

def get_geo(Ni, di, f_fit=False):
    """
    Calculate geometric mean diameter (GMD), geometric standard deviation (GSD), and optional fitting parameters.

    Parameters:
    Ni : array-like
        Particle number distribution.
    di : array-like
        Particle diameter distribution.
    f_fit : int, optional
        Flag for lognormal fitting. Defaults to 0 (use statistical definition).

    Returns:
    gmd : float
        Geometric mean diameter.
    gsd : float
        Geometric standard deviation.
    J : array-like, optional
        Output from the fitting optimization, only returned if f_fit is non-zero.
    """
    # Convert inputs to numpy arrays
    Ni = np.asarray(Ni, dtype=float)
    di = np.asarray(di, dtype=float)
    
    # Remove NaN values
    valid_indices = ~np.isnan(Ni) & ~np.isnan(di)
    Ni = Ni[valid_indices]
    di = di[valid_indices]

    # Calculate geometric mean diameter (GMD)
    gmd = np.exp(np.sum(Ni / np.nansum(Ni) * np.log(di)))
    gmd = np.real(gmd)

    # Calculate geometric standard deviation (GSD)
    log_gmd = np.log(gmd)
    gsd = np.exp(np.sqrt(np.sum(Ni / np.nansum(Ni) * (np.log(di) - log_gmd) ** 2)))
    gsd = np.real(gsd)

    # Cap lower limit of GSD, make monodisperse.
    if gsd < np.exp(0.06):
        gsd = 1

    if f_fit:
        # Lognormal fitting
        x0 = [np.max(Ni) * np.sqrt(2 * np.pi) * np.log(gsd), gmd, gsd]

        def fun(x):
            return Ni - norm.pdf(np.log(di), np.log(x[1]), np.log(x[2])) * np.sqrt(2 * np.pi) * np.log(x[2]) * x[0]

        result = least_squares(fun, x0, verbose=0)
        x1 = result.x
        J = result.cost
        gmd = x1[1]
        gsd = x1[2]
    else:
        J = None

    return gmd, gsd, J


def hc(mu, sg, q, a=None):
    """
    Compute new distribution moments using Hatch-Choate.

    Parameters:
    mu : float or array-like
        Geometric mean diameter (GMD) of the distribution.
    sg : float or array-like
        Geometric standard deviation (GSD) of the distribution.
    q : float or array-like
        Power used in the Hatch-Choate calculation (e.g., mass-mobility exponent).
    a : float or None, optional
        Pre-factor for the integrated variant. If None, computes the transformed mean.

    Returns:
    d : float or array-like
        The computed distribution moments.
    """
    mu = np.asarray(mu, dtype=float)
    sg = np.asarray(sg, dtype=float)
    q = np.asarray(q, dtype=float)

    if a is None:  # Transform the mean.
        d = mu * np.exp(q * (np.log(sg) ** 2))
    else:  # Integrated variant if pre-factor, a, is specified.
        a = np.asarray(a, dtype=float)
        d = a * mu ** q * np.exp(0.5 * q ** 2 * (np.log(sg) ** 2))

    return d


def dn2m(Ni, di, prop, type='HC'):
    """
    Compute PM concentration from a number distribution.

    Parameters:
    Ni : float or array-like
        Number concentration array.
    di : float or array-like
        Mobility diameter array.
    prop : dict or float
        Dictionary containing mass-mobility information.
    type : string or None, optional
        Type of conversion to perform:
            'HC' uses Hatch-Choate
            'HCF' uses Hatch-Chaote with lognormal fitting
            'NI' uses numerical integration.

    Returns:
    M : float or array-like
        Mass concentration.
    """
    di = di * 1e9

    if type == 'HC':
        dg, sg, _ = get_geo(Ni, di, False)
        H = hc(dg, sg, prop['zet'], prop['k'])
        M = H * np.nansum(Ni)  # output mass

    elif type == 'HCF':
        dg, sg, _ = get_geo(Ni, di, True)
        H = hc(dg, sg, prop['zet'], prop['k'])
        M = H * np.nansum(Ni)  # output mass

    elif type == 'NI':
        k = np.pi / 6 * prop['rho100'] * (100e-9 ** 3)
        M = k * np.nansum((di / 100) ** prop['zet'] * Ni)  # output mass
    
    return M


#== Other utilities =============================================#
def prop_update_flow(prop, Q):
    """
    Update the flow rate and v_bar in the properties dictionary.

    Parameters:
    prop (dict): Dictionary containing the properties.
    Q (float): Flow rate.

    Returns:
    dict: Updated properties dictionary.
    """
    prop['Q'] = Q
    prop['v_bar'] = prop['Q'] / prop['A']
    return prop

