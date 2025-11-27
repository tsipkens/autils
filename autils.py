import inspect
from collections import UserDict

import numpy as np
from scipy.optimize import fsolve, minimize, least_squares
from scipy.stats import norm

from autils import props

from tabulate import tabulate # used to show lists of dictionaries (e.g., for setpoint info)


def textdone():
    print('\r' +'\033[32m' + '^ DONE!' + '\033[0m' + '\n')


import numpy as np
from collections import UserDict
import inspect
from tabulate import tabulate  # Assuming this library is available for __repr__

class ComputedProperties(UserDict):
    """
    Generic container for properties. Subclasses collections.UserDict.
    
    - Data is stored in the self.data attribute (inherited from UserDict).
    - Supports vectorized inputs and dictionary-like/attribute access.
    """
    
    def __init__(self, **kwargs):
        # UserDict.__init__ initializes self.data = {}
        super().__init__()
        
        # Store properties directly in self.data (UserDict's internal dict)
        # This replaces the need for self._store
        self.data.update(kwargs)
        
        # 'toignore' must be stored as a regular attribute, not in self.data, 
        # so it doesn't get treated as a property.
        self.toignore = ['toignore']  
        
        # NOTE: We remove the manual __getitem__, __setitem__, __delitem__, 
        # __len__, __iter__, and __contains__ methods, as they are now 
        # correctly inherited and implemented by UserDict, operating on self.data.

    # --- Attribute Access (Simplified) ---
    # We must keep this custom __setattr__ to ensure new attributes are 
    # stored as properties in self.data, not as normal instance attributes.
    def __setattr__(self, name, value):
        # Allow setting of internal UserDict attributes and specific instance attributes
        if name in ('data', 'toignore'):
            super().__setattr__(name, value)
        # Treat all other assignments as property updates (i.e., put them in self.data)
        else:
            self.data[name] = value

    # We must keep this custom __getattr__ to allow attribute access (e.g., props.rc)
    # to retrieve items from the dictionary (self.data).
    def __getattr__(self, name):
        if name in self.data:
            return self.data[name]
        # This will raise AttributeError if the key isn't found in self.data
        raise AttributeError(f"'ComputedProperties' object has no attribute '{name}'")

    # --- Original Methods (Modified to use self.data) ---

    def check_lengths(self):
        """Broadcast all non-None attributes to the same length."""
        # Collect all non-None attributes except in toignore
        # Note: self.data is the dictionary store
        not_none = {k: np.atleast_1d(v) for k, v in self.data.items() 
                    if v is not None and k not in self.toignore}
        if not not_none:
            return
            
        target_len = max(v.size for v in not_none.values())
        for k, v in not_none.items():
            if v.size != target_len:
                self.data[k] = np.broadcast_to(v, target_len)
    
    def __repr__(self):
        """Display using tabulate."""
        keys = [k for k, v in self.data.items() if v is not None and k not in self.toignore]
        if not keys:
            return "<Empty Store>"
        
        if len(self.data[keys[0]]) == 1:
            """Override __repr__ to use self.data (UserDict internal store)."""
            lines = ["\r\033[32mProperties:\033[0m"]
            
            # Iterate over self.data, which holds all properties
            for attr, val in self.data.items():
                # Ensure 'toignore' is handled if it's stored in self.data for some reason
                if attr not in self.toignore:
                    lines.append(f"  \033[34m{attr}\033[0m → {repr(val)}")
            return "\n".join(lines)

        else:  # ensure all columns are the same size before stacking to form table
            try:
                cols = np.vstack([np.atleast_1d(self.data[k]).ravel() for k in keys]).T
            except ValueError:
                return "<Error: Property lengths do not match>"

            return tabulate(cols, headers=keys)

    def apply_functions(self, funcs):
        """Iteratively apply dictionary of functions to fill values."""
        values = {key: self.data.get(key, None) for key in funcs.keys()}
        changed = True
        
        while changed:
            changed = False
            for key, f in funcs.items():
                if values.get(key) is not None:
                    continue
                args = inspect.signature(f).parameters.keys()
                
                # Use self.data to check if required inputs are available
                if all(values.get(arg) is not None for arg in args):
                    # Compute the new value
                    computed_args = [values[arg] for arg in args]
                    values[key] = f(*computed_args)
                    changed = True
        
        # Reassign updated values (merges into self.data)
        self.data.update(values)
    
    # --- utilities for unique combinations ---
    def unique(self, tol=1e-4):
        """
        Return the unique combinations of property values across attributes.

        Parameters
        ----------
        tol : float
            Tolerance for considering values equal.

        Returns
        -------
        ComputedProperties
            Instance containing only unique combinations.
        np.ndarray
            Index mapping to reconstruct original instance.
        """
        keys = [k for k, v in self.data.items() if v is not None and not k in self.toignore]
        if not keys:
            return ComputedProperties(), np.array([], dtype=int)

        cols = np.vstack([np.atleast_1d(self.data[k]).ravel() for k in keys]).T

        # Compute magnitude and quantize for tolerance
        mags = np.floor(np.log10(np.abs(cols)))
        factor = 10**(-np.log10(tol) - mags - 1)
        quantized = np.round(cols * factor) / factor

        # Find unique rows and index mapping
        unique_rows, idx = np.unique(quantized, axis=0, return_inverse=True)
        unique_values = {k: unique_rows[:, ii] for ii, k in enumerate(keys)}

        return ComputedProperties(**unique_values), idx

    def expand(self, idx):
        """
        Reconstruct the original ComputedProperties instance
        from a unique instance and its index mapping.

        Parameters
        ----------
        idx : np.ndarray
            Indices of the unique rows in the original data.

        Returns
        -------
        ComputedProperties
            Reconstructed instance matching the original size.
        """
        keys = [k for k, v in self.data.items() if v is not None and not k in self.toignore]
        restored = {k: np.atleast_1d(self.data[k])[idx] for k in keys}
        return ComputedProperties(**restored)
    

class MassMob(ComputedProperties):
    def __init__(self, name=None, prop=None, **kwargs):
        """
        Initialize a MassMob object with material properties.

        Parameters
        ----------
        name : str, optional
            Name of a preset material. Default properties are loaded first.
        
        prop : dict, optional
            Dictionary of additional properties to store. Defaults to None.

        **kwargs : dict
            Arbitrary keyword arguments specifying properties. These
            override any preset values.
        """

        # --- Presets ---
        # Define presets dictionary. Use a mutable default value (dict) only for presets.
        presets = {
            "NaCl": {"zet": 3, "rho100": 2160},
            "salt": {"zet": 3, "rho100": 2160},
            "universal": {"zet": 2.48, "rho100": 510, "rhom": 1860, "dp100": 17.8, "DTEM": 0.35},
            "soot": {"zet": 2.48, "rho100": 510, "rhom": 1860, "dp100": 17.8, "DTEM": 0.35},
            "water": {"zet": 3, "rho100": 1000},
            "santovac": {"zet": 3, "rho100": 1198},
        }

        # 1. Prepare Initial Properties
        initial_props = {}
        if name and name in presets:
            # Load defaults from preset
            initial_props.update(presets[name])
        
        # 2. Add properties from the 'prop' dictionary (passed explicitly)
        if prop:
            initial_props.update(prop)
            
        # 3. Add properties from **kwargs (Highest precedence)
        initial_props.update(kwargs)

        # 4. Initialize Base Class (ComputedProperties)
        # This calls the parent's __init__, which populates self.data.
        super().__init__(**initial_props)

        # 5. Enforce Consistency and Compute Derived Properties
        self._solve()

    # ---- Enforce Consistency ----
    def _solve(self):
        """
        Defines and applies functions iteratively to fill and validate properties.
        Uses self.apply_functions inherited from ComputedProperties.
        """
        funcs = {
            'Dm': lambda zet: zet, 
            'zet': lambda Dm: Dm,
            # Note: 100e-9 is 100 nanometers
            'm100': lambda rho100: rho100 * np.pi / 6 * (100e-9) ** 3, 
            'm0': lambda m100, zet: m100 * (1 / 100) ** zet,
            'rho0': lambda m0: m0 * 6 / np.pi * 1e27,  # Conversion to rho0 (density at 1 nm, usually)
            'rho100': lambda m100: m100 * 6 / np.pi / (100e-9) ** 3,
            'k': lambda m0: m0  # alias
        }
        
        # apply_functions updates self.data directly, as implemented in the parent class
        self.apply_functions(funcs)


#== Functions for mass-mobility relations ======================#
def massmob_init(*args):
    """
    Fill in mass-mobility information using name-value pairs.
    Includes computing prop['m0'], which is used for mp2zp and mp2dm.
    (Partially deprecated. Not bridges to MassMob class for backwards compatibility.)

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
    
    d = {}  # initialize, then loop through arguments and format
    if isinstance(args[0], dict):
        d = args[0]
        _, *rest = args
        args = tuple(rest)
        d.pop('name')  # remove 'name', as will cause a conflict

    if len(args) > 1:
        for ii in range(0, len(args), 2):
            d[args[ii]] = args[ii+1]
    elif len(args) == 1:
        d['name'] = args[0]
    
    return MassMob(**d)


def massmob_add(prop, f1, v1=None, f2=None, v2=None):
    """
    Add or update mass-mobility parameters in an existing prop structure.

    Parameters:
    prop (dict or MassMob): The properties structure to update.
    f1 (str): Name of particle type OR field to add/update.
    v1 (any): Value corresponding to f1. If f1 is a particle type, this can be None.
    f2 (str): Second field to add/update (if using name-value pairs).
    v2 (any): Value corresponding to f2.

    Returns:
    MassMob: New MassMob instance with combined properties.
    """
    
    # 1. Extract existing properties from prop
    if isinstance(prop, MassMob):
        # Access the underlying dictionary from the MassMob instance
        existing_props = prop.data.copy()
    else:
        # If it's a regular dict, copy it
        existing_props = prop.copy()

    # 2. Determine new mass-mobility inputs
    if v1 is None:  # Case 1: f1 is a preset name (e.g., 'NaCl')
        inputs = {'name': f1}
    else:  # Case 2: f1, v1, f2, v2 are name-value pairs
        inputs = {f1: v1}
        if f2 is not None:
             inputs[f2] = v2
    
    prop = {**existing_props, **MassMob(**inputs)}  # merge properties
    prop = ComputedProperties(**prop)  # convert back to ComputedProperties

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


def newton(fun, x0, n=3):
    """
    Solve for a zero using Newtons method.
    """
    n = 3 if n is None else n  # set default
    h = x0 * 1e-5
    for _ in range(n):
        f = fun(x0)
        h = 1e-12
        fprime = (fun(x0 + h) - fun(x0)) / h
        x0 = x0 - f / fprime
    return x0

def fzero(fun, x0, n=3):
    """
    Wrapper to decide which zero solver to use.
    """
    if n == 0:
        return x0  # do not iterate
    elif n == np.inf:
        return fsolve(fun, x0)  # then use imported solver to default tolerances
    else:
        return newton(fun, x0, n)  # else do fast solve, with limited number of computations


#== Slip correction ============================================#
def cc(d, T=None, p=None, opt=None):
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
    if isinstance(T, str):
        opt = T
    elif T is not None:
        T = T

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

def mu(T=None, p=None):
    """
    Return gas viscosity [Pa*s], with two possible calculation methods.
    """
    if T is None or p is None:
        return 1.82e-5  # gas viscosity
    
    else:  # constants for Olfert laboratory / Kim et al.
        S = 110.4  # temperature [K]
        T_0 = 296.15  # reference temperature [K]
        vis_23 = 1.83245e-5  # reference viscosity [kg/(m*s)]
        return vis_23 * ((T / T_0) ** 1.5) * ((T_0 + S) / (T + S))  # gas viscosity

#== Aerodynamic diameter conversions ===========================#
def da2dm(da, prop, n=3, *args):
    """
    Convert aerodynamic diameter to mobility diameter.
    
    Parameters:
        da (float or numpy array): Aerodynamic diameter in meters.
        prop (dict): A dictionary containing the particle properties (e.g., density).
        n (int, optional): Number of iterations (see fzero).
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
    dm = fzero(fun_simple, da * 1e9, np.inf) * 1e-9

    # Iterative method (more precise)
    if n > 0:
        def fun_iter(dm):
            return 1e9 * (dm * 1e-9 * np.sqrt(dm2rhoeff(dm * 1e-9, prop) / rho0 *
                           cc(dm * 1e-9, *args) / cc(da, *args)) - da)
        dm = fzero(fun_iter, dm * 1e9, n) * 1e-9

    # Ensure result is real
    dm = np.real(dm)

    return dm

def da_rhoeff2dm(da, rho, n=3, *args):
    """
    Convert aerodynamic diameter and effective density to mobility diameter.
    
    Parameters:
        da (float or numpy array): Aerodynamic diameter in meters.
        rho (float or numpy array): Effective density of same size as da or scalar.
        n (int, optional): Number of iterations (also see fzero).
        *args: Additional arguments passed to the slip correction function Cc if needed.

    Returns:
        dm (float or numpy array): Mobility diameter in meters.
    """

    # Density of water in kg/m^3
    rho0 = 1e3

    # Calculate the mobility diameter directly without iteration (simple method)
    dm = da * np.sqrt(rho0 / rho)

    # Iterative method (more precise)
    if n > 0:
        def fun_iter(dm):
            return 1e9 * (dm * 1e-9 * np.sqrt(rho / rho0 *
                           cc(dm * 1e-9, *args) / cc(da, *args)) - da)
        
        dm = fzero(fun_iter, dm * 1e9, n) * 1e-9

    # Ensure result is real
    dm = np.real(dm)

    return dm

def da2dve(da, prop, n=3):
    """
    Convert the aerodynamic diameter to a volume-equivalent diameter.
    
    Parameters:
    da (float or ndarray): Aerodynamic diameter.
    prop (dict): Dictionary containing the material density and other properties.
    n (int, optional): Number of iterations (see also fzero).
    
    Returns:
    dve (float or ndarray): Volume-equivalent diameter.
    """
    rho0 = 1e3  # Density of water
    
    if n > 0:
        # Iterative method
        def fun(dve):
            return (dve * np.sqrt(prop['rhom'] / rho0 / dve2chi(dve, prop, n) * cc(dve) / cc(da)) - da) * 1e9
    else:
        # Simple method
        def fun(dve):
            return (da / np.sqrt(prop['rhom'] / rho0 / dve2chi(dve, prop, n)) - dve) * 1e9
    
    # Solve for dve iterativekly.
    dve = fzero(fun, da, n)
    
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
def dm2chi(dm, prop, n=3, *args):
    """
    Compute the dynamic shape factor at a given mobility diameter.

    Parameters:
    dm: Mobility diameter
    prop: Dictionary with properties including 'rhom'
    n (int, optional): Number of iterations (also see fzero).
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
    if n > 0:
        chi = chi / cc(dm, *args) * cc(dve, *args)

    # Alternative form (commented in MATLAB code)
    # b = (6 / np.pi * prop['k'] * 1e9 ** prop['zet'] / prop['rhom']) ** (1 / 3)
    # chi = dm ** (1 - prop['zet'] / 3) / b
    # if n:
    #     chi = chi / cc(dm) * cc(b * dm ** (prop['zet'] / 3))

    return chi


def dm2da(dm, prop, n=3, *args):
    """
    Convert the mobility diameter (dm) to an aerodynamic diameter (da).

    Parameters:
    dm (float or ndarray): Mobility diameter.
    prop (dict): Dictionary containing properties for conversion.
    n (int, optional): Number of iterations (also see fzero).
    args: Additional arguments for the Cunningham slip correction factor (Cc).

    Returns:
    da (float or ndarray): Aerodynamic diameter.
    """
    rho0 = 1e3  # density of water in kg/m^3

    # Compute effective density
    rhoeff = dm2rhoeff(dm, prop)  # effective density from dm
    da = dm * np.sqrt(rhoeff / rho0)  # simple aerodynamic diameter

    # Alternate, iterative method
    if n > 0:
        def func(da1):
            return 1e9 * (da * np.sqrt(cc(dm, *args) / cc(da1, *args)) - da1)

        # Solve for aerodynamic diameter by finding zero.
        da = fzero(func, da, n)

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


def dm2mp(dm, prop):
    """
    Computes the effective density from mass and mobility diameter.
    
    Parameters:
        dm: float or array-like
            Mobility diameter in meters.
        prop: dict
            Properties containing relevant physical parameters.
    
    Returns:
        m: float or array-like
            Particle mass in kilograms.
    """
    return prop['m0'] * (dm * 1e9) ** prop['zet']


def dm2rhoeff(dm, prop):
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
    return rhoeff(dm, dm2mp(dm, prop))


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
    zp (float or ndarray): Electromobility.
    """
    # Define constants
    e = 1.6022e-19  # electron charge [C]

    # Default gas viscosity for Buckley/Davies
    B = cc(dm,T,p) / (3 * np.pi * mu(T,p) * dm)  # mechanical mobility

    # Calculate electromobility
    zp = B * e * z

    return B, zp


def dm2dpp(dm, prop):
    """
    Computed primary particle size from a da-dpp relationship, where da = dm.
    """
    return (prop['dp100'] * 1e-9) * (dm / 100e-9) ** prop['DTEM']


def dm2npp(dm, prop):
    """
    Computed mobility diameter from the number of primary particles and a combination of
    (1) the mass-mobility relation and (2) a da-dpp relationship, where da = dm.

    NOTE: Requires that prop contains both the standard mass-mobility parameters as
    well as a relation for the primary particle size. 
    """

    N100 = prop['rho100'] / prop['rhom'] * (100 / prop['dp100']) ** 3
    Npp = N100 * (dm / 100e-9) ** (3 * (1 - prop['DTEM']))

    return Npp


def npp2dm(npp, prop):
    """
    Computed number of primary particles from mobility diameter and a combination of
    (1) the mass-mobility relation and (2) a da-dpp relationship, where da = dm.

    NOTE: Requires that prop contains both the standard mass-mobility parameters as
    well as a relation for the primary particle size. 
    """

    N100 = prop['rho100'] / prop['rhom'] * (100 / prop['dp100']) ** 3

    return 100e-9 * (npp / N100) ** (1 / (3 * (1 - prop['DTEM'])))


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


def dm_rhoeff2da(dm, rho, n=3, *args):
    """
    Convert aerodynamic diameter and effective density to mobility diameter.
    
    Parameters:
        da (float or numpy array): Aerodynamic diameter in meters.
        rho (float or numpy array): Effective density of same size as da or scalar.
        n (int, optional): Number of iterations (also see fzero).
        *args: Additional arguments passed to the slip correction function Cc if needed.

    Returns:
        dm (float or numpy array): Mobility diameter in meters.
    """

    # Density of water in kg/m^3
    rho0 = 1e3

    # Calculate the mobility diameter directly without iteration (simple method)
    da = dm * np.sqrt(rho / rho0)

    # Iterative method (more precise)
    if n > 0:
        def fun_iter(da):
            return 1e9 * (dm * np.sqrt(rho / rho0 *
                           cc(dm, *args) / cc(da * 1e-9, *args)) - da * 1e-9)
        
        da = fzero(fun_iter, da * 1e9, n) * 1e-9

    # Ensure result is real
    da = np.real(da)

    return da


def dm_mp2da(dm, mp, n=3, *varargs):
    """
    Compute aerodynamic diameter from mobility diameter and mass.

    Parameters:
    dm (float or ndarray): Mobility diameter in meters.
    mp (float or ndarray): Particle mass in kilograms.
    n (int, optional): Number of iterations (also see fzero).
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    da (float or ndarray): Aerodynamic diameter in meters.
    """

    # Compute simple volume-equivalent and aerodynamic diameters.
    da = dm**(-1/2) * mp**(1/2) * np.sqrt(6 / (np.pi * 1e3))

    if n > 0:
        # Define the function to solve
        def fun(da1):
            return 1e9 * (da * np.sqrt(cc(dm, *varargs) / cc(da1, *varargs)) - da1)

        # Solve for aerodynamic diameter
        da = fzero(fun, da, n)

    return da


#== Volume-equivalent diameter conversions =======================#
def dve2chi(dve, prop, n=3):
    """
    Compute the dynamic shape factor at a given volume-equivalent diameter (dve).
    
    Parameters:
    dve (float or ndarray): Volume-equivalent diameter.
    prop (dict): Dictionary containing properties like 'rhom'.
    n (int, optional): Number of iterations (also see fzero).
    
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
    if n > 0:
        chi = (dm / dve) * (cc(dve) / cc(dm))

    return chi


def dve2da(dve, prop, n=3, *varargs):
    """
    Convert volume-equivalent diameter to aerodynamic diameter.

    Parameters:
    dve (float or ndarray): Volume-equivalent diameter in meters.
    prop (dict): Dictionary containing properties including 'rhom'.
    n (int, optional): Number of iterations (also see fzero).
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    da (float or ndarray): Aerodynamic diameter in meters.
    """
    rho0 = 1e3  # Density of water in kg/m^3

    # Compute simple volume-equivalent and aerodynamic diameters
    chi = dve2chi(dve, prop, n)
    da = dve * np.sqrt(prop['rhom'] / rho0 / chi)

    if n > 0:
        # Define the function to solve
        def fun(da1):
            return 1e9 * (dve * np.sqrt(prop['rhom'] / rho0 / chi * 
                                         cc(dve, *varargs) / cc(da1, *varargs)) - da1)

        # Solve for aerodynamic diameter
        da = fzero(fun, da, n)

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
def mp2dm(mp, prop):
    """
    Calculate mobility diameter from particle mass using mass-mobility relation.

    Parameters:
    mp (float or ndarray): Particle mass in kilograms.
    prop (dict): Dictionary containing mass-mobility properties with 'm0' and 'zet'.

    Returns:
    float or ndarray: Mobility diameter in meters.
    """
    if 'm0' not in prop:
        prop = massmob_init(prop)

    # Compute mobility diameter and return
    return 1e-9 * (mp / prop['m0']) ** (1 / prop['zet'])


def mp_da2dm(mp, da, n=3, *varargs):
    """
    Compute mobility diameter from particle mass and aerodynamic diameter.

    Parameters:
    mp (float or ndarray): Particle mass in kilograms.
    da (float or ndarray): Aerodynamic diameter in meters.
    n (int, optional): Number of iterations (also see fzero).
    *varargs: Optional arguments for Cunningham slip correction calculation.

    Returns:
    dm (float or ndarray): Mobility diameter in meters.
    """
    # Compute simple volume-equivalent and aerodynamic diameters
    dm = da**(-2) * mp * (6 / (np.pi * 1e3))

    if n > 0:
        # Define the function to solve
        def fun(dm1):
            return 1e9 * (dm * (cc(dm1, *varargs) / cc(da, *varargs)) - dm1)

        # Solve for mobility diameter
        dm = fzero(fun, dm, n)

    return dm


def mp2zp(mp, prop, z=1, T=None, p=None):
    """
    Calculate electric mobility from particle mass.

    Parameters:
    mp (float or ndarray): Particle mass in kilograms.
    z (int or ndarray): Integer charge state.
    T (float, optional): Temperature in Kelvin.
    p (float, optional): Pressure in atm.
    prop (dict): Dictionary containing mass-mobility properties with 'm0' and 'zet'.

    Returns:
    tuple: Mechanical mobility (B), electromobility (zp), and mobility diameter (d).
    """
    # Compute mobility diameter
    d = mp2dm(mp, prop)

    # Compute mechanical and electrical mobility
    B, zp = dm2zp(d, z, T, p)

    return B, zp, d


#== Mobility > mobility diameter ================================#
def zp2dm(zp, z=1, T=None, p=None):
    """
    Calculate mobility diameter from a vector of mobilities.
    
    Parameters:
    zp : array-like
        The electric mobility values.
    z : int
        Integer charge state.
    T : float, optional
        Temperature in Kelvin. If not specified, defaults to using Buckley/Davies method.
    p : float, optional
        Pressure in atm. If not specified, defaults to using Buckley/Davies method.
        
    Returns:
    dm : array-like
        The mobility diameters.
    B : array-like
        Mechanical mobility.
    """
    e = 1.6022e-19  # Electron charge [C]
    B = zp / (e * z)  # Mechanical mobility
    
    dm0 = 1 / (3 * np.pi * mu(T,p) * B)  # Initial estimate
    dm0 = dm0 * cc(dm0,T,p)

    def fun(dm):
        return B - cc(np.abs(dm),T,p) / (3 * np.pi * mu(T,p) * np.abs(dm))
    
    dm0 = fzero(fun, dm0, n=100)
    
    return dm0, B


#== Distribution manipulation ===================================#
def sdm2sda(sd, d, prop, n=3):
    """
    Calculate the GSD for the aerodynamic distribution from mobility distribution GSD.

    Parameters:
    sd (float or ndarray): GSD for the mobility distribution.
    d (float or ndarray): Mobility diameter.
    prop (dict): Dictionary containing mass-mobility properties.
    n (int, optional): Number of iterations (also see fzero).

    Returns:
    float or ndarray: GSD for the aerodynamic distribution.
    """
    d1 = np.exp(np.log(d) * 1.01)
    d2 = np.exp(np.log(d) * 0.99)

    da1 = dm2da(d1, prop, n)
    da2 = dm2da(d2, prop, n)

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

