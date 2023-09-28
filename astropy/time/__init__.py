# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.time`.
    """

    use_fast_parser = _config.ConfigItem(
        ["True", "False", "force"],
        "Use fast C parser for supported time strings formats, including ISO, "
        "ISOT, and YearDayTime. Allowed values are the 'False' (use Python parser),"
        "'True' (use C parser and fall through to Python parser if fails), and "
        "'force' (use C parser and raise exception if it fails). Note that the"
        "options are all strings.",
    )

    masked_array_type = _config.ConfigItem(
        ["astropy", "numpy"],
        'The type of masked array used for masked output data.  Can be "astropy" '
        'for `astropy.utils.masked.Masked` or "numpy" to use `numpy.ma.MaskedArray`. '
        "Note that if `astropy.units.Quantity` is produced, the output always "
        "uses `astropy.utils.masked.Masked`, since `numpy.ma.MaskedArray` does not "
        "work with quantities.",
    )
    # Create a dict of available masked classes for speed.
    # Use local imports so we do not pollute the module namespace.
    from numpy.ma import MaskedArray

    from astropy.utils.masked import Masked

    _MASKED_CLASSES = {"astropy": Masked, "numpy": MaskedArray}

    @property
    def _masked_cls(self):
        """The masked class set by ``masked_array_type``.

        This is |Masked| for "astropy", `numpy.ma.MaskedArray` for "numpy".
        """
        return self._MASKED_CLASSES[self.masked_array_type]


conf = Conf()

# isort: off
from .formats import *
from .core import *

# isort: on
