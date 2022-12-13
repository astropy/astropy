# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.table`.
    """

    use_fast_parser = _config.ConfigItem(
        ["True", "False", "force"],
        "Use fast C parser for supported time strings formats, including ISO, "
        "ISOT, and YearDayTime. Allowed values are the 'False' (use Python parser),"
        "'True' (use C parser and fall through to Python parser if fails), and "
        "'force' (use C parser and raise exception if it fails). Note that the"
        "options are all strings.",
    )


conf = Conf()

# isort: off
from .formats import *
from .core import *

# isort: on
