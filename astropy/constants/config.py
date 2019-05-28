# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Configures the codata and iaudata used, possibly using user configuration.
"""
# Note: doing this in __init__ causes import problems with units,
# as si.py and cgs.py have to import the result.
from astropy import physical_constants, astronomical_constants

if ((physical_constants.get() == 'codata2018') or
        (physical_constants.get() == 'astropyconst40')):
    from . import codata2018 as codata
elif ((physical_constants.get() == 'codata2014') or
        (physical_constants.get() == 'astropyconst20')):
    from .astropyconst20 import codata2014 as codata  # noqa
elif ((physical_constants.get() == 'codata2010') or
        (physical_constants.get() == 'astropyconst13')):
    from .astropyconst13 import codata2010 as codata  # noqa
else:
    # ScienceState validates values so this should never happen
    raise ValueError('Invalid physical constants version: {}'
                     .format(physical_constants.get()))

if ((astronomical_constants.get() == 'iau2015') or
        (astronomical_constants.get() == 'astropyconst40') or
        (astronomical_constants.get() == 'astropyconst20')):
    from . import iau2015 as iaudata
elif ((astronomical_constants.get() == 'iau2012') or
        (astronomical_constants.get() == 'astropyconst13')):
    from .astropyconst13 import iau2012 as iaudata  # noqa
else:
    # ScienceState validates values so this should never happen
    raise ValueError('Invalid astronomical constants version: {}'
                     .format(astronomical_constants.get()))
