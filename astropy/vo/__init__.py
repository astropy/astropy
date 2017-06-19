# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The `~astropy.vo` subpackage provides virtual observatory (VO) related
functionality.
"""

from .. import config as _config


# NOTE: This is deprecated along with other Cone Search stuff, but it feels
#       weird for config item to be issuing deprecation warnings.

class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.vo`.
    """

    vos_baseurl = _config.ConfigItem(
        'http://stsdas.stsci.edu/astrolib/vo_databases/',
        'URL where the VO Service database file is stored.',
        aliases=['astropy.vo.client.vos_catalog.vos_baseurl'])
    conesearch_dbname = _config.ConfigItem(
        'conesearch_good',
        'Conesearch database name to use.',
        aliases=['astropy.vo.client.conesearch.conesearch_dbname'])


conf = Conf()
