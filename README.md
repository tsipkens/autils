# AUTILS

(***A***erosol ***Util***itie***s***)

A set of utilities for computing and comparing aerosol properties. This code is still under construction and could change at any time. 

### General notes

Unless otherwise stated, quantities are given in SI units. So:

+ Diameters are given in m. 
+ Masses are given in kg.
+ Densities are given in kg/m<sup>3</sup>.  

### +massmob package

This package contains a set of tools for generating and manipulating mass-mobility property structures. The mass-mobility exponent is stored in the `zet` field, but a duplicate value is stored in a `Dm` field for backwards compatibility. This prop structure can then be given to functions like `dm2mp(...)` to make use of the underlying mass-mobility relation parameters. 
