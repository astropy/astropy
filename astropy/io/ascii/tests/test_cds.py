# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some methods related to ``CDS`` format
reader/writer.
Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""
from io import StringIO

from astropy.table import Table


def test_write_data():
    output = ['S08-229 4625 1.23 1.23 -1.5      ',
              'S05-10  4342 0.91 1.82 -2.11 0.14',
              'S05-47  4654 1.28 1.74 -1.64 0.16']

    dat = 'data/cdsFunctional2.dat'
    t = Table.read(dat, format='ascii.cds')

    out = StringIO()
    t.write(out, format='ascii.cds')
    assert out.getvalue().splitlines() == output
