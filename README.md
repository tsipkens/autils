# AUTILS

(***A***erosol ***Util***itie***s***)

A set of utilities for computing and comparing aerosol properties. 

For tutorial files and UQ tools, see the **tutorial** branch. The main branch has been streamlined for addition as a submodule in other code bases. 

A JavaScript/web-enabled version of some of this code is available at [tsipkens/autils-web](https://github.com/tsipkens/autils-web), which was deployed as a web app at [tsipkens.github.io/autils-web](https://tsipkens.github.io/autils-web/). 

## General notes

Unless otherwise stated, quantities are given in SI units. So:

+ Diameters are given in m. 
+ Masses are given in kg.
+ Densities are given in kg/m<sup>3</sup>.

Many of the conversions also depend on a **prop** structure, a MATLAB struct containing fields that describe properties of the material. Most functions only require the mass-mobility relationship, which can be initiated using the functions in the **+massmob** package, described below. In addition to these properties, one can also add: 

+ The material desnity, `prop.rhom`. 
+ The quantities `prop.Dtem` and `prop.dp100`, which are used to compute a primary particle diameter from a mobility diameter. 
+ The fractal prefactor, `prop.kf`, and fractal dimension, `prop.Df`, which are used to compute a radius of gyration. 

## Packages

### +massmob

This package contains a set of tools for generating and manipulating mass-mobility property structures. The mass-mobility exponent is stored in the `zet` field, but a duplicate value is stored in a `Dm` field for backwards compatibility. This prop structure can then be given to functions like `dm2mp(...)` to make use of the underlying mass-mobility relation parameters. The structure, when intiated using `massmob.init(...)`, will contain all of the renditions of the parameters related to the mass-mobility power law. 

### +uq

This package can be used for error modeling, including computing covariances and simulating artificial signals.

### +pm

A package with tools for computing total particle mass (PM) concentration from number-based measurements, e.g., using numerical integration.
