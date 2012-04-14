# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains dictionaries with sets of parameters for a
given cosmology.

Each cosmology has the following parameters defined:

    Oc          Omega cold dark matter
    Ob          Omega baryon
    Om          Omega matter
    Ol          Omega lambda
    H0          Hubble parameter at z=0 in km/s/Mpc
    n           Density perturbation spectral index
    sigma8      Density perturbation amplitude
    tau         Ionisation optical depth
    z_reion     Redshift of hydrogen reionisation
    t0          Age of the universe in Gyr

The list of cosmologies available are given by the tuple
`available`. Current cosmologies available:

WMAP 7 year parameters from Komatsu et al. 2011, ApJS, 192, 18. Table
1 (WMAP + BAO + H0 ML).

WMAP 5 year parameters from Komatsu et al. 2009, ApJS, 180, 330. Table
1 (WMAP + BAO + SN ML).

Both these cosmologies are flat (omega matter + omega lambda = 1).
"""

# Komatsu et al. 2011, WMAP + BAO + H0 ML (table 1).

WMAP7 = dict(
    Oc = 0.226,
    Ob = 0.0455,
    Om = 0.272,
    Ol = 0.728,
    H0 = 70.4,
    n = 0.967,
    sigma8 = 0.810,
    tau = 0.085,
    z_reion = 10.3,
    t0 = 13.76,
    )

# Komatsu et al. 2009 WMAP + BAO + SN ML (table 1).

WMAP5 = dict(
    Oc = 0.231,
    Ob = 0.0459,
    Om = 0.277,
    Ol = 0.723,
    H0 = 70.2,
    n = 0.962,
    sigma8 = 0.817,
    tau = 0.088,
    z_reion = 11.3,
    t0 = 13.72,
    )

available = tuple(k for k in locals() if not k.startswith('_'))
