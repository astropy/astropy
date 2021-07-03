# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some methods related to ``CDS`` format
reader/writer.
Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""
from io import StringIO

from astropy.table import Table
from astropy.io import ascii
from astropy import units as u


def test_write_data():
    output = ['S08-229 4625 1.23 1.23 -1.50      ',
              'S05-10  4342 0.91 1.82 -2.11  0.14',
              'S05-47  4654 1.28 1.74 -1.64  0.16']

    dat = 'data/cdsFunctional2.dat'
    t = Table.read(dat, format='ascii.cds')

    out = StringIO()
    t.write(out, format='ascii.cds')
    # get the last section of table which will be the data.
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines)
             if s.startswith(('------', '======='))]
    lines = lines[i_secs[-1]+1:]
    assert lines == output


test_dat = ['names e d s i',
            'HD81809 1E-7 22.25608 +2 67',
            'HD103095 -31.6e5 +27.2500 -9E34 -30']


def test_write_ByteByByte_units():
    t = ascii.read(test_dat)
    colUnits = [None, u.C, u.kg, None, u.year]
    t._set_column_attribute('unit', colUnits)
    # add a column with magnitude units.
    # note that magnitude has to be assigned for each value explicitly.
    t['magnitude'] = [u.Magnitude(25), u.Magnitude(-9)]
    colUnits.append(u.mag)
    out = StringIO()
    t.write(out, format='ascii.cds')
    # read written table.
    tRead = ascii.read(out.getvalue(), format='cds')
    assert [tRead[col].unit for col in tRead.columns] == colUnits
