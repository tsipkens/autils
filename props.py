
import os

import numpy as np

import re
import yaml


#== General property constructors ==========================#
def covert_prop_array(prop):
    """
    Coverts all entries in prop to lists. 
    This is required later for flexible transfer function evaluation.
    """
    for key, value in prop.items():
        if not type(prop[key]) == list:
            prop[key] = np.array([value])

    return prop


def pma(opts='cpma'):
    prop = {}

    if opts == 'olfert':
        opts = 'cpma'  # for legacy code

    with open(os.path.join(os.path.dirname(__file__), 'config/pma', opts + '.yaml')) as file:
        # Explicitly specify loader.
        # Allows for more variability in scientific notation.
        loader = yaml.SafeLoader
        loader.add_implicit_resolver(
            u'tag:yaml.org,2002:float',
            re.compile(u'''^(?:
             [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
            |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
            |\\.[0-9_]+(?:[eE][-+][0-9]+)?
            |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
            |[-+]?\\.(?:inf|Inf|INF)
            |\\.(?:nan|NaN|NAN))$''', re.X),
            list(u'-+0123456789.'))

        prop = yaml.load(file, Loader=loader)

    if prop == {}:
        print('Specified PMA property set is not available. Reverting to default.')
        prop = {
            'r1': 0.06,  # inner electrode radius [m]
            'r2': 0.061,  # outer electrode radius [m]
            'L': 0.2,  # length of chamber [m]
            'p': 1,  # pressure [atm]
            'T': 293,  # system temperature [K]
            'Q': 3 / 1000 / 60,  # volume flow rate (m ** 3/s) (prev: ~1 lpm)
            'omega_hat': 32 / 33  # ratio of angular speeds
        }

    # -- Parameters related to CPMA geometry ----------------------------------#
    prop['rc'] = (prop['r1'] + prop['r2']) / 2
    prop['r_hat'] = prop['r1'] / prop['r2']
    prop['del'] = (prop['r2'] - prop['r1']) / 2  # half gap width

    prop['A'] = np.pi * (prop['r2'] ** 2 - prop['r1'] ** 2)  # cross sectional area of APM
    prop['v_bar'] = prop['Q'] / prop['A']  # average flow velocity

    # -- For diffusion --------------------------------------------------------#
    kB = 1.3806488e-23  # Boltzmann's constant
    prop['D'] = lambda B: kB * prop['T'] * B  # diffusion coefficient

    # -- Default mass-mobility information -------------#
    prop['Dm'] = 3  # mass-mobility exponent
    prop['zet'] = prop['Dm']
    prop['m0'] = 4.7124e-25  # mass-mobility relation density
    # Common alternate: Dm = 2.48 m0 = 2.9280e-24

    return covert_prop_array(prop)


def update_flow(prop, Q):
    """
    PROP_UPDATE_FLOW  A simple utility to update the flow rate and v_bar.
    Acts as a reminder that these two parameters are linked.
    To be used in conjunction with particle mass analyzer prop structures.

    AUTHOR: Timothy Sipkens, 2021-03-26
    """
    prop['Q'] = Q

    if 'v_bar' in prop.keys():
        prop['v_bar'] = prop['Q'] / prop['A']

    return prop
    

def dma(*args):
    # Parse inputs
    if len(args) == 0:
        opts = {'param': ''}
    elif len(args) == 1:
        if isinstance(args[0], dict):
            opts = args[0]
        else:
            opts = {'param': args[0]}
    elif len(args) == 2:
        opts = {'param': '', 'model': args[0], 'flow': args[1]}
    else:
        raise ValueError('Too many inputs.')

    # Parameter specs. Default is first option.
    opts['param'] = opts.get('param', '')
    if opts['param'] == 'olfert':
        opts['model'] = '3018'
        opts['flow'] = 'low'
    elif opts['param'] == 'buckley':
        opts['model'] = '3081'
        opts['flow'] = 'low'
    else:
        if opts['param']:
            opts['model'] = opts['param']  # assume param is just the model number

    # Assign defaults
    opts['model'] = opts.get('model', '3018')
    opts['flow'] = opts.get('flow', 'low')

    prop = {
        'model': opts['model'],
        'flow': opts['flow'],
        'T': 293,  # default temperature [K]
        'p': 1,    # default pressure
    }

    # Assign classifier dimensions
    if prop['model'] in ['3018', '3081', '3080']:
        prop['L'] = 0.44369  # length, m
        prop['R1'] = 0.00937  # inner radius, m
        prop['R2'] = 0.01961  # outer radius, m
    elif prop['model'] == 'custom':
        prop.update(opts.get('prop', {}))  # read in from opts

    # Assign flow rates
    if prop['flow'] == 'low':
        prop['Qs'] = 0.3 / 60 / 1000  # Sample flow [m^3/s]
        prop['Qa'] = 0.3 / 60 / 1000  # Aerosol flow [m^3/s]
        prop['Qc'] = 3 / 60 / 1000    # Sheath flow [m^3/s]
        prop['Qm'] = 3 / 60 / 1000    # Exhaust flow [m^3/s]
    elif prop['flow'] == 'high':
        prop['Qs'] = 1.5 / 60 / 1000
        prop['Qa'] = 1.5 / 60 / 1000
        prop['Qc'] = 15 / 60 / 1000
        prop['Qm'] = 15 / 60 / 1000
    elif prop['flow'] == 'buckley':
        prop['Qc'] = 4.89E-3 / 60  # sheath flow [m^3/s]
        prop['Qa'] = 1.02E-3 / 60  # aerosol flow [m^3/s]
        prop['Qs'] = prop['Qa']    # sample flow [m^3/s]
        prop['Qm'] = prop['Qc']    # exhaust flow [m^3/s]

    return covert_prop_array(prop)


def aac():
    # Initialize the property dictionary
    prop = {}

    # Default temperature and pressure
    prop['T'] = 293  # default temperature [K]
    prop['p'] = 1    # default pressure [atm]

    # Classifier dimensions
    prop['r1'] = 0.056  # 0.043;
    prop['r2'] = 0.060  # 0.045;
    prop['L'] = 0.206   # 0.21;

    # Low flow
    prop['Qs'] = 0.3 / 60 / 1000  # sample flow (leaving classifier) [m^3/s]
    prop['Qa'] = prop['Qs']       # aerosol flow (entering classifier) [m^3/s]
    prop['Qsh'] = 3 / 60 / 1000   # sheath flow [m^3/s]
    prop['Qexh'] = prop['Qsh'] + prop['Qa'] - prop['Qs']  # outlet sheath [m^3/s]

    # Uncomment the following lines for high flow:
    # prop['Qs'] = 1.5 / 60 / 1000  # sample flow (leaving classifier) [m^3/s]
    # prop['Qa'] = prop['Qs']       # aerosol flow (entering classifier) [m^3/s]
    # prop['Qsh'] = 15 / 60 / 1000  # sheath flow [m^3/s]
    # prop['Qexh'] = prop['Qsh'] + prop['Qa'] - prop['Qs']  # outlet sheath [m^3/s]

    return covert_prop_array(prop)


def elpi(spec=None):
    """
    Get properties of impactors, such as the ELPI.
    
    Parameters:
    spec (str): Specification type. Default is 'default'.
    
    Returns:
    dict: A dictionary containing d50 and s50 values.
    """
    if spec is None:
        spec = 'default'

    prop = {}

    if spec in ['default', 'elpi']:
        # Values from Järvinen et al. (2014).
        prop['d50'] = [15.7, 30.4, 54.1, 94.3, 154., 
                       254., 380., 600., 943., 1620, 
                       2460, 3640, 5340]
        prop['s50'] = [3.32, 3.65, 3.89, 3.05, 3.62, 
                       6.30, 8.43, 7.16, 6.21, 5.32, 
                       5.33, 4.14, 3.66]

    return prop
