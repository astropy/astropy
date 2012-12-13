# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Temporary solution until `astropy.vo.server.cs_mstr_list`
includes `testQuery` fields.

"""
from __future__ import print_function, division

# STDLIB
from xml.dom import minidom

# LOCAL
from ...utils import OrderedDict  # For 2.6 compatibility
from ...utils.data import get_readable_fileobj


def parse_cs(id):
    """Return `testQuery` pars as dict for given Resource ID."""
    url = 'http://nvo.stsci.edu/vor10/getRecord.aspx?' \
        'id={}&format=xml'.format(id)
    tqp = ['ra', 'dec', 'sr']
    d = OrderedDict()
    with get_readable_fileobj(url) as fd:
        dom = minidom.parse(fd)
    tq = dom.getElementsByTagName('testQuery')
    for key in tqp:
       d[key.upper()] = \
           tq[0].getElementsByTagName(key)[0].firstChild.nodeValue.strip()
    return d
