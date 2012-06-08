"""
Tools to verify the syntax of unit attributes according to the
`Standards for Astronomical Catalogues, Version 2.0
<http://cdsarc.u-strasbg.fr/doc/catstd-3.2.htx>`_
"""

# STDLIB
import re

# LOCAL
from .exceptions import warn_or_raise, W50


__all__ = ['unit_names', 'unit_prefixes', 'is_unit', 'check_unit']


unit_names = u"""
    - % a A AU arcmin arcsec barn bit byte C cd ct D d deg eV F g h H
    Hz J Jy K lm lx m mag mas min mol N Ohm Pa pc pix rad Ry s S
    solLum solMass solRad Sun sr T V W Wb yr
    """.split()


unit_prefixes = u"""
    d c m u n p f a z y da h k M G T P E Z Y
    """.split()


_unit_name_regex = u'|'.join(u'(?:{0})'.format(x) for x in unit_names)
_unit_prefix_regex = u'|'.join(u'(?:{0})'.format(x) for x in unit_prefixes)


def is_unit(s):
    """
    Returns `True` if *s* is a valid unit string as defined by `Standards for Astronomical Catalogues, Version 2.0 <http://cdsarc.u-strasbg.fr/doc/catstd-3.2.htx>`_.
    """
    number = ur'[+\-]?[0-9]+(?:\.[0-9]*)?(?:[+\-][0-9]+)?'
    factor = ur'{0}(?:x{0})?'.format(number)
    unit = ur'(?:{0})?(?:{1})'.format(_unit_prefix_regex, _unit_name_regex)
    bracketed_unit = ur'(?:{0})|(?:\[{0}\])'.format(unit)
    expression = (
        ur'^(?P<factor>{0})?(?P<unit>{1})(?P<subunit>[/.]{1})?(?P<power>{2})?$'.format(
        factor, bracketed_unit, number))

    match = re.match(expression, s)
    return match is not None


def check_unit(unit, attr_name, config={}, pos=None):
    """
    Raises a `ValueError` if *unit* is not a valid unit as defined by `Standards for Astronomical Catalogues, Version 2.0 <http://cdsarc.u-strasbg.fr/doc/catstd-3.2.htx>`_.
    """
    if (unit is not None and not is_unit(unit)):
        warn_or_raise(W50, W50, unit, config, pos)
        return False
    return True
