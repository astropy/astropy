# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Temporary solution until `astropy.vo.validator.validate.CS_MSTR_LIST`
includes ``<testQuery>`` fields.

"""
from __future__ import print_function, division

# STDLIB
from xml.dom import minidom

# LOCAL
from ...utils import OrderedDict  # For 2.6 compatibility
from ...utils.data import get_readable_fileobj


def parse_cs(id):
    """Return ``<testQuery>`` pars as dict for given Resource ID."""
    if isinstance(id, bytes):  # pragma: py3
        id = id.decode('ascii')
    url = 'http://nvo.stsci.edu/vor10/getRecord.aspx?' \
        'id={0}&format=xml'.format(id)
    tqp = ['ra', 'dec', 'sr']
    d = OrderedDict()
    with get_readable_fileobj(url, encoding='binary') as fd:
        dom = minidom.parse(fd)
    tq = dom.getElementsByTagName('testQuery')
    for key in tqp:
        d[key.upper()] = \
            tq[0].getElementsByTagName(key)[0].firstChild.nodeValue.strip()
    return d
