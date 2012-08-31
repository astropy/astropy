# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities used by a number of the formats.
"""


def get_grouped_by_powers(unit):
    """
    Groups the powers and bases in the given
    `~astropy.units.core.CompositeUnit` into positive powers and
    negative powers for easy display on either side of a solidus.

    Parameters
    ----------
    unit : `astropy.units.core.CompositeUnit`

    Returns
    -------
    positives, negatives : tuple of lists
       Each element in each list is tuple of the form (*base*,
       *power*).  The negatives have the sign of their power reversed
       (i.e. the powers are all positive).
    """
    positive = []
    negative = []
    for base, power in zip(unit.bases, unit.powers):
        if power < 0:
            negative.append((base, -power))
        elif power > 0:
            positive.append((base, power))
        else:
            raise ValueError("Unit with 0 power")
    return positive, negative


def split_mantissa_exponent(v):
    """
    Given a number, split it into its mantissa and base 10 exponent
    parts, each as strings.
    """
    parts = []

    x = "{0:.2e}".format(val).split('e')
    if x[0] != '1.00':
        m = x[0]
    else:
        m = ''

    ex = x[1].lstrip("0+")
    if len(ex) > 0 and ex[0] == '-':
        ex = '-' + ex[1:].lstrip('0')

    return m, ex
