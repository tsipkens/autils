# AUTILS

(***A***erosol ***Util***itie***s***)

A set of utilities for computing and comparing aerosol properties. This code is still under construction and could change at any time.

A JavaScript/web-enabled version of some of this code is available at [tsipkens/autils-web](https://github.com/tsipkens/autils-web), which was deployed as a web app at [tsipkens.github.io/autils-web](https://tsipkens.github.io/autils-web/). 

## General notes

Unless otherwise stated, quantities are given in SI units. So:

+ Diameters are given in m. 
+ Masses are given in kg.
+ Densities are given in kg/m<sup>3</sup>.

## Packages

### +massmob

This package contains a set of tools for generating and manipulating mass-mobility property structures. The mass-mobility exponent is stored in the `zet` field, but a duplicate value is stored in a `Dm` field for backwards compatibility. This prop structure can then be given to functions like `dm2mp(...)` to make use of the underlying mass-mobility relation parameters.

### +uq

This package can be used for error modeling, including computing covariances and simulating artificial signals.

### +pm

A package with tools for computing total particle mass (PM) concentration from number-based measurements, e.g., using numerical integration.
