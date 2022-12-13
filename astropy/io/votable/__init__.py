# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This package reads and writes data formats used by the Virtual
Observatory (VO) initiative, particularly the VOTable XML format.
"""


from astropy import config as _config

from .exceptions import (
    IOWarning,
    UnimplementedWarning,
    VOTableChangeWarning,
    VOTableSpecError,
    VOTableSpecWarning,
    VOWarning,
)
from .table import from_table, is_votable, parse, parse_single_table, validate, writeto

__all__ = [
    "Conf",
    "conf",
    "parse",
    "parse_single_table",
    "validate",
    "from_table",
    "is_votable",
    "writeto",
    "VOWarning",
    "VOTableChangeWarning",
    "VOTableSpecWarning",
    "UnimplementedWarning",
    "IOWarning",
    "VOTableSpecError",
]


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.io.votable`.
    """

    verify = _config.ConfigItem(
        "ignore",
        "Can be 'exception' (treat fixable violations of the VOTable spec as "
        "exceptions), 'warn' (show warnings for VOTable spec violations), or "
        "'ignore' (silently ignore VOTable spec violations)",
        aliases=["astropy.io.votable.table.pedantic", "astropy.io.votable.pedantic"],
    )


conf = Conf()
