
import inspect
import re
import numpy as np
from autils import autils

dm = np.array([150e-9])

# Test two calls to MassMob.
prop = autils.MassMob(zet=3, rho100=510)
prop = autils.MassMob(name='soot')

# Initial load.
params = {
    'z': 1, 'dm': dm, 
    'mp': autils.dm2mp(dm, prop),
    'rho': autils.dm2rhoeff(dm, prop),
    'da': autils.dm2da(dm, prop),
    'dve': autils.dm2dve(dm, prop),
    'chi': autils.dm2chi(dm, prop),
    'dpp': autils.dm2dpp(dm, prop),
    'npp': autils.dm2npp(dm, prop)}
params['zp'] = autils.dm2zp(dm, z=params['z'])[1]
params['prop'] = prop
params['n'] = 3

# Initiate structure of outputs with initial values from above.
outputs = params.copy()
for key in ['prop', 'z', 'n']:
    outputs.pop(key)
outputs['rhoeff'] = outputs['rho']
outputs.pop('rho')


for name, func in inspect.getmembers(autils, inspect.isfunction):
    if "2" not in name:
        continue

    print(f"\nTesting: {name}")

    # Example: parse "mp2zp" -> inputs ["mp"], outputs ["zp"]
    match = re.match(r"([a-z]+)2([a-z]+)", name)
    if match:
        inp, out = match.groups()
        inps = np.array([inp])

    else:
        match = re.match(r"([a-z]+)_([a-z]+)2([a-z]+)", name)
        if match:
            inp1, inp2, out = match.groups()
            inps = np.array([inp1, inp2])

        else:
            print(f"  Skipping (name does not fit pattern): {name}")
            continue

    inps[inps == 'rhoeff'] = 'rho'  # catch this exception
    print(f"  Expected input: {inps}, output: {out}")

    if not all([inp in params.keys() for inp in inps]):
        print(f"  \033[93mSkipping (some inputs are not defined): {name}\033[0m")
        continue

    # # Inspect function signature for args
    sig = inspect.signature(func)
    plist = list(sig.parameters.keys())

    # Get arguments.
    args = []
    for p in plist:
        if p in params.keys():
            args.append(params[p])
        else:
            break
        
    # Print full list of inputs, include optional arguments.
    print(f"  Full input: {plist}")

    # Try calling function and catch errors. 
    try:
        result = func(*args)
        print(f"  \033[92mSuccess. Output: {out}={result}\033[0m")
    except Exception as e:
        print(f"  \033[91mError calling {name}: {e}\033[0m")
        next
    
    try:
        if isinstance(result, tuple):
            outputs[out] = np.append(outputs[out], result[1])
        else:
            outputs[out] = np.append(outputs[out], result)
    except:
        print(f"  Warning: Output not saved.")


print('\nOUTPUTS=')
for key in outputs.keys():
    print(f"{key}={outputs[key]}")