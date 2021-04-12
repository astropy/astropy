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

- Planck 2018 (Planck18) parameters from Planck Collaboration 2020,
  A&A, 641, A6 (Paper VI), Table 2 (TT, TE, EE + lowE + lensing + BAO)

- Planck 2018 (Planck18_arXiv_v2) parameters from Planck Collaboration 2018,
  arXiv:1807.06209v2 (Paper VI), Table 2 (TT, TE, EE + lowE + lensing + BAO)

- Planck 2015 (Planck15) parameters from Planck Collaboration 2016, A&A, 594,
  A13 (Paper XIII), Table 4 (TT, TE, EE + lowP + lensing + ext)

- Planck 2013 (Planck13) parameters from Planck Collaboration 2014, A&A, 571,
  A16 (Paper XVI), Table 5 (Planck + WP + highL + BAO)

- WMAP 9 year (WMAP9) parameters from Hinshaw et al. 2013, ApJS, 208, 19,
  doi: 10.1088/0067-0049/208/2/19. Table 4 (WMAP9 + eCMB + BAO + H0)

- WMAP 7 year (WMAP7) parameters from Komatsu et al. 2011, ApJS, 192, 18,
  doi: 10.1088/0067-0049/192/2/18. Table 1 (WMAP + BAO + H0 ML).

- WMAP 5 year (WMAP5) parameters from Komatsu et al. 2009, ApJS, 180, 330,
  doi: 10.1088/0067-0049/180/2/330. Table 1 (WMAP + BAO + SN ML).

Note: if you add a new cosmology, please also update the table
in the 'Built-in Cosmologies' section of astropy/docs/cosmology/index.rst
in addition to the list above.  You also need to add them to the ``available``
variable.

"""

import json
import pathlib
import sys
from types import MappingProxyType

__all__ = ["available"]

# Read in the built-in realizations' parameter dicts file
with open(pathlib.Path(__file__).parent.joinpath("parameters.json"), 'r') as f:
    _parameter_dict = json.loads(f.read())

# register all realizations' parameter dicts as "available"
available = frozenset(_parameter_dict.keys())

# Add each parameter set to this module as an immutable dictionary.
for name, params in _parameter_dict.items():
    # set parameter dictionary here
    setattr(sys.modules[__name__], name, MappingProxyType(params))
    __all__.append(name)

# clear namespace
del name, params
