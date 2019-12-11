# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This subpackage contains classes and functions for celestial coordinates
of astronomical objects. It also contains a framework for conversions
between coordinate systems.
"""

from .errors import *
from .angles import *
from .baseframe import *
from .attributes import *
from .distances import *
from .earth import *
from .transformations import *
from .builtin_frames import *
from .name_resolve import *
from .matching import *
from .representation import *
from .sky_coordinate import *
from .funcs import *
from .calculation import *
from .solar_system import *


from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.coordinates`.
    """

    skycoord_init_counter_warn_threshold = _config.ConfigItem(
        100,
        'This controls the (globally-counted) number of SkyCoord '
        'initializations that will trigger a warning to the user. To disable '
        'this counter and warning, set this config item to 0.'
    )


conf = Conf()
