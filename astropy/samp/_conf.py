# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.samp`.
    """

    use_internet = _config.ConfigItem(
        True,
        "Whether to allow `astropy.samp` to use the internet, if available.",
        aliases=["astropy.samp.utils.use_internet"],
    )

    n_retries = _config.ConfigItem(
        10, "How many times to retry communications when they fail"
    )
conf = Conf()
