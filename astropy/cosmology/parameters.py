# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains dictionaries with sets of parameters for a
given cosmology.

Each cosmology has the following parameters defined:

    ==========  =====================================
    Oc          Omega cold dark matter
    Ob          Omega baryon
    Om          Omega matter
    Ode         Omega dark energy
    w0          Equation of state of dark energy at z=0 (Pressure/density)
    wa          a derivative of equation of state
    H0          Hubble parameter at z=0 in km/s/Mpc
    n           Density perturbation spectral index
    sigma8      Density perturbation amplitude
    tau         Ionisation optical depth
    z_reion     Redshift of hydrogen reionisation
    t0          Age of the universe in Gyr
    reference   Reference for the parameters
    ==========  =====================================

The list of cosmologies available are given by the tuple
`available`. Current cosmologies available:

WMAP 7 year parameters from Komatsu et al. 2011, ApJS, 192, 18. Table
1 (WMAP + BAO + H0 ML).

WMAP 5 year parameters from Komatsu et al. 2009, ApJS, 180, 330. Table
1 (WMAP + BAO + SN ML).

Both these cosmologies are flat (omega matter + omega dark energy = 1)
and assume dark enegy is the cosmological constant (w0 = -1, wa = 0).
"""

# Komatsu et al. 2011, WMAP + BAO + H0 ML (table 1).

WMAP7 = dict(
    Oc  = 0.226,
    Ob  = 0.0455,
    Om  = 0.272,
    Ode = 0.728,
    w0  = -1.0,
    wa  = 0.0,
    H0 = 70.4,
    n = 0.967,
    sigma8 = 0.810,
    tau = 0.085,
    z_reion = 10.3,
    t0 = 13.76,
    reference = ("Komatsu et al. 2011, ApJS, 192, 18. "
                 "Table 1 (WMAP + BAO + H0 ML)")
    )

# Komatsu et al. 2009 WMAP + BAO + SN ML (table 1).

WMAP5 = dict(
    Oc  = 0.231,
    Ob  = 0.0459,
    Om  = 0.277,
    Ode = 0.723,
    w0  = -1.0,
    wa  = 0.0,
    H0 = 70.2,
    n = 0.962,
    sigma8 = 0.817,
    tau = 0.088,
    z_reion = 11.3,
    t0 = 13.72,
    reference = ("Komatsu et al. 2009, ApJS, 180, 330. "
                 "Table 1 (WMAP + BAO + SN ML)")
    )

available = tuple(k for k in locals() if not k.startswith('_'))
