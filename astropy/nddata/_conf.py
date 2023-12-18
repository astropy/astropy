# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.nddata`.
    """

    warn_unsupported_correlated = _config.ConfigItem(
        True,
        "Whether to issue a warning if `~astropy.nddata.NDData` arithmetic "
        "is performed with uncertainties and the uncertainties do not "
        "support the propagation of correlated uncertainties.",
    )

    warn_setting_unit_directly = _config.ConfigItem(
        True,
        "Whether to issue a warning when the `~astropy.nddata.NDData` unit "
        "attribute is changed from a non-``None`` value to another value "
        "that data values/uncertainties are not scaled with the unit change.",
    )
conf = Conf()
