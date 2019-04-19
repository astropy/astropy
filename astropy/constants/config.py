# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Configures the codata and iaudata used, possibly using user configuration.
"""
# Note: doing this in __init__ causes import problems with units,
# as si.py and cgs.py have to import the result.
from . import conf

if ((conf.physical_constants == 'codata2014') or
        (conf.physical_constants == 'astropyconst20')):
    from . import codata2014 as codata
elif ((conf.physical_constants == 'codata2010') or
        (conf.physical_constants == 'astropyconst13')):
    from .astropyconst13 import codata2010 as codata  # noqa
else:
    raise ValueError('Invalid physical constants version: {}'
                     .format(conf.physical_constants))

if ((conf.astronomical_constants == 'iau2015') or
        (conf.astronomical_constants == 'astropyconst20')):
    from . import iau2015 as iaudata
elif ((conf.astronomical_constants == 'iau2012') or
        (conf.astronomical_constants == 'astropyconst13')):
    from .astropyconst13 import iau2012 as iaudata  # noqa
else:
    raise ValueError('Invalid astronomical constants version: {}'
                     .format(conf.astronomical_constants))
