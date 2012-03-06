""" Cosmological parameters from the 5-year and 7-year WMAP results.

WMAP 7 year parameters are from Komatsu et al. 2011, ApJS, 192,
18. Table 1 (WMAP + BAO + H0 ML).

WMAP 5 year parameters are from Komatsu et al. 2009, ApJS, 180,
330. Table 1 (WMAP + BAO + SN ML).

Each set of parameters is in a dictionary with the following keys:

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

Both of these cosmologies are flat (omega matter + omega lambda = 1).
"""

# Komatsu et al. 2011, WMAP + BAO + H0 ML (table 1).
#
# This is a flat cosmology (Om + Ol = 1)

WMAP7 = dict(
    Oc =      0.226,  # Omega cold dark matter
    Ob =      0.0455, # Omega baryon
    Om =      0.272,  # Omega matter
    Ol =      0.728,  # Omega lambda
    H0 =      70.4,   # Hubble parameter at z=0
    n =       0.967,  # spectral index
    sigma8 =  0.810,  # sigma8
    tau =     0.085,  # 
    z_reion = 10.3,   # Redshift of hydrogen reionisation
    t0 =      13.76,  # Age of the universe
    )

# Komatsu et al. 2009 WMAP + BAO + SN ML (table 1). Also flat.

WMAP5 = dict(
    Oc =      0.231,  # Omega cold dark matter
    Ob =      0.0459, # Omega baryon
    Om =      0.277,  # Omega matter
    Ol =      0.723,  # Omega lambda
    H0 =      70.2,   # Hubble parameter at z=0
    n =       0.962,  # spectral index
    sigma8 =  0.817,  # sigma8
    tau =     0.088,  # 
    z_reion = 11.3,   # Redshift of hydrogen reionisation
    t0 =      13.72,  # Age of the universe
    )
