# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Configures the codata and iaudata used, possibly using user configuration.
"""
# Note: doing this in __init__ causes import problems with units,
# as si.py and cgs.py have to import the result.
import importlib

import astropy

phys_version = astropy.physical_constants.get()
astro_version = astropy.astronomical_constants.get()

codata = importlib.import_module(".constants." + phys_version, "astropy")
iaudata = importlib.import_module(".constants." + astro_version, "astropy")
