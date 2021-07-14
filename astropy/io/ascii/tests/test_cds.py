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

# flake8: noqa: E291

test_dat = ['names e d s i',
            'HD81809 1E-7 22.25608 +2 67',
            'HD103095 -31.6e5 +27.2500 -9E34 -30']


def test_write_data():
    exp_output = ['S08-229 4625 1.23 1.23 -1.50      ',
                  'S05-10  4342 0.91 1.82 -2.11  0.14',
                  'S05-47  4654 1.28 1.74 -1.64  0.16']

    dat = 'data/cdsFunctional2.dat'
    t = Table.read(dat, format='ascii.cds')

    out = StringIO()
    t.write(out, format='ascii.cds')
    # Get the last section of table which will be the data.
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines)
              if s.startswith(('------', '======='))]
    lines = lines[i_secs[-1]+1:]
    assert lines == exp_output


def test_write_byte_by_byte_units():
    t = ascii.read(test_dat)
    col_units = [None, u.C, u.kg, None, u.year]
    t._set_column_attribute('unit', col_units)
    # Add a column with magnitude units.
    # Note that magnitude has to be assigned for each value explicitly.
    t['magnitude'] = [u.Magnitude(25), u.Magnitude(-9)]
    col_units.append(u.mag)
    out = StringIO()
    t.write(out, format='ascii.cds')
    # Read written table.
    tRead = ascii.read(out.getvalue(), format='cds')
    assert [tRead[col].unit for col in tRead.columns] == col_units


def test_write_readme_with_default_options():
    exp_output = '''\
                                                      (1st author ?, Date ?)
================================================================================
Title ?
    Authors ?
    ref ?
================================================================================
Keywords: 

Objects:
    -----------------------------------------
       RA   (2000)   DE    Designation(s)
    -----------------------------------------

Abstract:
  Abstract ?

File Summary:
--------------------------------------------------------------------------------
 FileName    Lrecl   Records    Explanations
--------------------------------------------------------------------------------
ReadMe         80        . this file
table          35        2          


--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 8   A8     ---    names    Description of names             
10-14   E5.1   ---    e       [-3160000.0/0.01] Description of e
16-23   F8.5   ---    d       [22.25/27.25] Description of d    
25-31   E7.1   ---    s       [-9e+34/2.0] Description of s     
33-35   I3     ---    i       [-30/67] Description of i         

--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by author  / pyreadme )
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67
HD103095 -3e+06 27.25000  -9e+34 -30
'''
    t = ascii.read(test_dat)
    out = StringIO()
    t.write(out, format='ascii.cds')
    assert out.getvalue() == exp_output


def test_write_empty_table():
    out = StringIO()
    import pytest
    with pytest.raises(ValueError):
        Table().write(out, format='ascii.cds')


def test_write_null_data_values():
    exp_output = ['HD81809  1e-07   22.25608  2.0e+00  67',
                  'HD103095 -3e+06  27.25000 -9.0e+34 -30',
                  'Sun                        5.3e+27    ']
    t = ascii.read(test_dat)
    t.add_row(['Sun', '3.25', '0', '5.3e27', '2'],
                mask=[False, True, True, False, True])
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines)
              if s.startswith(('------', '======='))]
    lines = lines[i_secs[-1]+1:]
    assert lines == exp_output


def test_write_byte_by_byte_for_masked_column():
    """
    This test differs from the ``test_write_null_data_values``
    above in that it tests the column value limits in the Byte-By-Byte
    description section for columns whose values are masked.
    It also checks the description for columns with same values.
    """
    exp_output = '''\
--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 8   A8     ---    names    Description of names         
10-14   E5.1   ---    e       [0.0/0.01]? Description of e  
16-18   F3.0   ---    d       ? Description of d            
20-26   E7.1   ---    s       [-9e+34/2.0] Description of s 
28-30   I3     ---    i       [-30/67] Description of i     
32-34   F3.1   ---    sameF   [5.0/5.0] Description of sameF
36-37   I2     ---    sameI   [20] Description of sameI     

--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by author  / pyreadme )
--------------------------------------------------------------------------------
HD81809  1e-07    2e+00  67 5.0 20
HD103095         -9e+34 -30 5.0 20
'''
    from astropy.table import MaskedColumn
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name='sameF')
    t.add_column([20, 20], name='sameI')
    t['e'] = MaskedColumn(t['e'], mask=[False, True])
    t['d'] = MaskedColumn(t['d'], mask=[True, True])
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines)
              if s.startswith(('------', '======='))]
    lines = lines[i_secs[-6]:]
    assert lines == exp_output.splitlines()
