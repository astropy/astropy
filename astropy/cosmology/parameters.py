# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" This module contains dictionaries with sets of parameters for a
given cosmology.

Each cosmology has the following parameters defined:

    ==========  =====================================
    Oc0         Omega cold dark matter at z=0
    Ob0         Omega baryon at z=0
    Om0         Omega matter at z=0
    flat        Is this assumed flat?  If not, Ode0 must be specified
    Ode0        Omega dark energy at z=0 if flat is False
    H0          Hubble parameter at z=0 in km/s/Mpc
    n           Density perturbation spectral index
    Tcmb0       Current temperature of the CMB
    Neff        Effective number of neutrino species
    sigma8      Density perturbation amplitude
    tau         Ionisation optical depth
    z_reion     Redshift of hydrogen reionisation
    t0          Age of the universe in Gyr
    reference   Reference for the parameters
    ==========  =====================================

The list of cosmologies available are given by the tuple
`available`. Current cosmologies available:

Planck 15 parameters from Planck Collaboration 2015, arXiv: 1502.01589
 (Paper XIII), Table 4 (TT, TE, EE + lowP + lensing + ext)

Planck13 parameters from Planck Collaboration 2013, arXiv:1303.5076
 (Paper XVI), Table 5 (Planck + WP + highL + BAO)

WMAP 9 year parameters from Hinshaw et al. 2013, ApJS, 208, 19,
doi: 10.1088/0067-0049/208/2/19. Table 4 (WMAP9 + eCMB + BAO + H0)

WMAP 7 year parameters from Komatsu et al. 2011, ApJS, 192, 18,
doi: 10.1088/0067-0049/192/2/18. Table 1 (WMAP + BAO + H0 ML).

WMAP 5 year parameters from Komatsu et al. 2009, ApJS, 180, 330,
doi: 10.1088/0067-0049/180/2/330. Table 1 (WMAP + BAO + SN ML).

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
# delete these things from the namespace so we can automatically find
# all of the parameter dictionaries below.
del absolute_import
del division
del print_function
del unicode_literals

# Note: if you add a new cosmology, please also update the table
# in the 'Built-in Cosmologies' section of astropy/docs/cosmology/index.rst
# in addition to the list above.  You also need to add them to
# __all__ in core.py

# Planck 2015 paper XII Table 4 final column (best fit)
Planck15 = dict(
    Oc0=0.2589,
    Ob0=0.04860,
    Om0=0.3075,
    H0=67.74,
    n=0.9667,
    sigma8=0.8159,
    tau=0.066,
    z_reion=8.8,
    t0=13.799,
    Tcmb0=2.7255,
    Neff=3.046,
    flat=True,
    m_nu=[0., 0., 0.06],
    reference=("Planck Collaboration 2015, Paper XII, arXiv:1502.01589"
               " Table 4 (TT, TE, EE + lowP + lensing + ext)")
)

# Planck 2013 paper XVI Table 5 penultimate column (best fit)
Planck13 = dict(
    Oc0=0.25886,
    Ob0=0.048252,
    Om0=0.30712,
    H0=67.77,
    n=0.9611,
    sigma8=0.8288,
    tau=0.0952,
    z_reion=11.52,
    t0=13.7965,
    Tcmb0=2.7255,
    Neff=3.046,
    flat=True,
    m_nu=[0., 0., 0.06],
    reference=("Planck Collaboration 2013, Paper XVI, arXiv:1303.5076"
               " Table 5 (Planck + WP + highL + BAO)")
)


WMAP9 = dict(
    Oc0=0.2402,
    Ob0=0.04628,
    Om0=0.2865,
    H0=69.32,
    n=0.9608,
    sigma8=0.820,
    tau=0.081,
    z_reion=10.1,
    t0=13.772,
    Tcmb0=2.725,
    Neff=3.04,
    m_nu=0.0,
    flat=True,
    reference=("Hinshaw et al. 2013, ApJS, 208, 19, "
               "doi: 10.1088/0067-0049/208/2/19. "
               "Table 4 (WMAP9 + eCMB + BAO + H0, last column)")
)

WMAP7 = dict(
    Oc0=0.226,
    Ob0=0.0455,
    Om0=0.272,
    H0=70.4,
    n=0.967,
    sigma8=0.810,
    tau=0.085,
    z_reion=10.3,
    t0=13.76,
    Tcmb0=2.725,
    Neff=3.04,
    m_nu=0.0,
    flat=True,
    reference=("Komatsu et al. 2011, ApJS, 192, 18, "
               "doi: 10.1088/0067-0049/192/2/18. "
               "Table 1 (WMAP + BAO + H0 ML).")
)

WMAP5 = dict(
    Oc0=0.231,
    Ob0=0.0459,
    Om0=0.277,
    H0=70.2,
    n=0.962,
    sigma8=0.817,
    tau=0.088,
    z_reion=11.3,
    t0=13.72,
    Tcmb0=2.725,
    Neff=3.04,
    m_nu=0.0,
    flat=True,
    reference=("Komatsu et al. 2009, ApJS, 180, 330, "
               "doi: 10.1088/0067-0049/180/2/330. "
               "Table 1 (WMAP + BAO + SN ML).")
)

available = tuple(k for k in locals() if not k.startswith('_'))
