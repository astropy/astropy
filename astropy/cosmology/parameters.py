# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains dictionaries with sets of parameters for a
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
    m_nu        Assumed mass of neutrino species, in eV.
    sigma8      Density perturbation amplitude
    tau         Ionisation optical depth
    z_reion     Redshift of hydrogen reionisation
    t0          Age of the universe in Gyr
    reference   Reference for the parameters
    ==========  =====================================

The list of cosmologies available are given by the tuple
`available`. Current cosmologies available:

Planck 2018 (Planck18) parameters from Planck Collaboration 2020,
 A&A, 641, A6 (Paper VI), Table 2 (TT, TE, EE + lowE + lensing + BAO)

Planck 2018 (Planck18_arXiv_v2) parameters from Planck Collaboration 2018,
 arXiv:1807.06209v2 (Paper VI), Table 2 (TT, TE, EE + lowE + lensing + BAO)

Planck 2015 (Planck15) parameters from Planck Collaboration 2016, A&A, 594, A13
 (Paper XIII), Table 4 (TT, TE, EE + lowP + lensing + ext)

Planck 2013 (Planck13) parameters from Planck Collaboration 2014, A&A, 571, A16
 (Paper XVI), Table 5 (Planck + WP + highL + BAO)

WMAP 9 year (WMAP9) parameters from Hinshaw et al. 2013, ApJS, 208, 19,
doi: 10.1088/0067-0049/208/2/19. Table 4 (WMAP9 + eCMB + BAO + H0)

WMAP 7 year (WMAP7) parameters from Komatsu et al. 2011, ApJS, 192, 18,
doi: 10.1088/0067-0049/192/2/18. Table 1 (WMAP + BAO + H0 ML).

WMAP 5 year (WMAP5) parameters from Komatsu et al. 2009, ApJS, 180, 330,
doi: 10.1088/0067-0049/180/2/330. Table 1 (WMAP + BAO + SN ML).

WMAP 3 year (WMAP3) parameters from Spergel et al. 2007, ApJS, 170, 377,
doi:  10.1086/513700. Table 6. (WMAP + SNGold) Obtained from https://lambda.gsfc.nasa.gov/product/map/dr2/params/lcdm_wmap_sngold.cfm
Tcmb0 and Neff are the standard values as also used for WMAP5, 7, 9.
Pending WMAP team approval and subject to change.

WMAP 1 year (WMAP1) parameters from Spergel et al. 2003, ApJS, 148, 175,
doi:  10.1086/377226. Table 7 (WMAP + CBI + ACBAR + 2dFGRS + Lya)
Tcmb0 and Neff are the standard values as also used for WMAP5, 7, 9.
Pending WMAP team approval and subject to change.

"""

# STDLIB
from types import MappingProxyType

# LOCAL
import astropy.units as u

# Note: if you add a new cosmology, please also update the table
# in the 'Built-in Cosmologies' section of astropy/docs/cosmology/index.rst
# in addition to the list above.  You also need to add them to the 'available'
# list at the bottom of this file.

# Planck 2018 paper VI
# Unlike Planck 2015, the paper includes massive neutrinos in Om0, which here
# are included in m_nu.  Hence, the Om0 value differs slightly from the paper.
Planck18 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.2607,
    Ob0=0.04897,
    Om0=0.30966,
    H0=67.66 * (u.km / u.s / u.Mpc),
    n=0.9665,
    sigma8=0.8102,
    tau=0.0561,
    z_reion=7.82,
    t0=13.787 * u.Gyr,
    Tcmb0=2.7255 * u.K,
    Neff=3.046,
    flat=True,
    m_nu=[0., 0., 0.06] * u.eV,
    reference=("Planck Collaboration 2018, 2020, A&A, 641, A6  (Paper VI),"
               " Table 2 (TT, TE, EE + lowE + lensing + BAO)")
))

# Planck 2018 paper VI v2.  Identical to Planck18 above.
# Warning: deprecated and will be removed in future versions.
Planck18_arXiv_v2 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.2607,
    Ob0=0.04897,
    Om0=0.30966,
    H0=67.66 * (u.km / u.s / u.Mpc),
    n=0.9665,
    sigma8=0.8102,
    tau=0.0561,
    z_reion=7.82,
    t0=13.787 * u.Gyr,
    Tcmb0=2.7255 * u.K,
    Neff=3.046,
    flat=True,
    m_nu=[0., 0., 0.06] * u.eV,
    reference=("DEPRECATED: Planck Collaboration 2018, arXiv:1807.06209 v2 (Paper VI),"
               " Table 2 (TT, TE, EE + lowE + lensing + BAO)")
))

# Planck 2015 paper XII Table 4 final column (best fit)
Planck15 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.2589,
    Ob0=0.04860,
    Om0=0.3075,
    H0=67.74 * (u.km / u.s / u.Mpc),
    n=0.9667,
    sigma8=0.8159,
    tau=0.066,
    z_reion=8.8,
    t0=13.799 * u.Gyr,
    Tcmb0=2.7255 * u.K,
    Neff=3.046,
    flat=True,
    m_nu=[0., 0., 0.06] * u.eV,
    reference=("Planck Collaboration 2016, A&A, 594, A13 (Paper XIII),"
               " Table 4 (TT, TE, EE + lowP + lensing + ext)")
))

# Planck 2013 paper XVI Table 5 penultimate column (best fit)
Planck13 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.25886,
    Ob0=0.048252,
    Om0=0.30712,
    H0=67.77 * (u.km / u.s / u.Mpc),
    n=0.9611,
    sigma8=0.8288,
    tau=0.0952,
    z_reion=11.52,
    t0=13.7965 * u.Gyr,
    Tcmb0=2.7255 * u.K,
    Neff=3.046,
    flat=True,
    m_nu=[0., 0., 0.06] * u.eV,
    reference=("Planck Collaboration 2014, A&A, 571, A16 (Paper XVI),"
               " Table 5 (Planck + WP + highL + BAO)")
))


WMAP9 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.2402,
    Ob0=0.04628,
    Om0=0.2865,
    H0=69.32 * (u.km / u.s / u.Mpc),
    n=0.9608,
    sigma8=0.820,
    tau=0.081,
    z_reion=10.1,
    t0=13.772 * u.Gyr,
    Tcmb0=2.725 * u.K,
    Neff=3.04,
    m_nu=0.0 * u.eV,
    flat=True,
    reference=("Hinshaw et al. 2013, ApJS, 208, 19, "
               "doi: 10.1088/0067-0049/208/2/19. "
               "Table 4 (WMAP9 + eCMB + BAO + H0, last column)")
))

WMAP7 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.226,
    Ob0=0.0455,
    Om0=0.272,
    H0=70.4 * (u.km / u.s / u.Mpc),
    n=0.967,
    sigma8=0.810,
    tau=0.085,
    z_reion=10.3,
    t0=13.76 * u.Gyr,
    Tcmb0=2.725 * u.K,
    Neff=3.04,
    m_nu=0.0 * u.eV,
    flat=True,
    reference=("Komatsu et al. 2011, ApJS, 192, 18, "
               "doi: 10.1088/0067-0049/192/2/18. "
               "Table 1 (WMAP + BAO + H0 ML).")
))

WMAP5 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.231,
    Ob0=0.0459,
    Om0=0.277,
    H0=70.2 * (u.km / u.s / u.Mpc),
    n=0.962,
    sigma8=0.817,
    tau=0.088,
    z_reion=11.3,
    t0=13.72 * u.Gyr,
    Tcmb0=2.725 * u.K,
    Neff=3.04,
    m_nu=0.0 * u.eV,
    flat=True,
    reference=("Komatsu et al. 2009, ApJS, 180, 330, "
               "doi: 10.1088/0067-0049/180/2/330. "
               "Table 1 (WMAP + BAO + SN ML).")
))

WMAP3 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.230,
    Ob0=0.0454,
    Om0=0.276,
    H0=70.1 * (u.km / u.s / u.Mpc),
    n=0.946,
    sigma8=0.784,
    tau=0.079,
    z_reion=10.3,
    t0=13.78 * u.Gyr,
    Tcmb0=2.725 * u.K,
    Neff=3.04,
    m_nu=0.0 * u.eV,
    flat=True,
    reference=(r"Spergel et al. 2007, ApJS, 170, 377, "
               r"doi:  10.1086/513700. "
               r"Table 6 (WMAP + SNGold) "
               r"obtained from: https://lambda.gsfc.nasa.gov/product/map/dr2/params/lcdm_wmap_sngold.cfm"
               r"\nPending WMAP team approval and subject to change.")
))

WMAP1 = MappingProxyType(dict(
    cosmology="FlatLambdaCDM",
    Oc0=0.213,
    Ob0=0.0436,
    Om0=0.257,
    H0=72. * (u.km / u.s / u.Mpc),
    n=0.96,
    sigma8=0.75,
    tau=0.117,
    z_reion=17.0,  # Only from WMAP1. Does not exist in the combined analysis.
    t0=13.4 * u.Gyr,
    Tcmb0=2.725 * u.K,
    Neff=3.04,
    m_nu=0.0 * u.eV,
    flat=True,
    reference=(r"Spergel et al. 2003, ApJS, 148, 175, "
               r"doi:  10.1086/377226. "
               r"Table 7 (WMAP + CBI + ACBAR + 2dFGRS + Lya)."
               r"\nPending WMAP team approval and subject to change.")
))


# If new parameters are added, this list must be updated
available = ('Planck18', 'Planck18_arXiv_v2', 'Planck15', 'Planck13',
             'WMAP9', 'WMAP7', 'WMAP5', 'WMAP3', 'WMAP1')
