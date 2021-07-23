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


test_dat = ['names e d s i',
            'HD81809 1E-7 22.25608 +2 67',
            'HD103095 -31.6e5 +27.2500 -9E34 -30']


def test_write_data():
    exp_output = ['S05-5   4337 0.77 1.80 -2.07      ',
                  'S08-229 4625 1.23 1.23 -1.50      ',
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
    col_units = [None, u.C, u.kg, u.m/u.s, u.year]
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
 1- 8  A8     ---    names   Description of names              
10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e
16-23  F8.5   ---    d       [22.25/27.25] Description of d    
25-31  E7.1   ---    s       [-9e+34/2.0] Description of s     
33-35  I3     ---    i       [-30/67] Description of i         
--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by 1st author ?  / astropy.io.ascii )
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67
HD103095 -3e+06 27.25000  -9e+34 -30
''' # noqa: W291
    t = ascii.read(test_dat)
    out = StringIO()
    t.write(out, format='ascii.cds')
    assert out.getvalue() == exp_output


def test_write_empty_table():
    out = StringIO()
    import pytest
    with pytest.raises(NotImplementedError):
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
    lines = lines[i_secs[-1]+1:]  # Last section is the data.
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
 1- 8  A8     ---    names   Description of names          
10-14  E5.1   ---    e       [0.0/0.01]? Description of e  
16-18  F3.0   ---    d       ? Description of d            
20-26  E7.1   ---    s       [-9e+34/2.0] Description of s 
28-30  I3     ---    i       [-30/67] Description of i     
32-34  F3.1   ---    sameF   [5.0/5.0] Description of sameF
36-37  I2     ---    sameI   [20] Description of sameI     
--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by 1st author ?  / astropy.io.ascii )
--------------------------------------------------------------------------------
HD81809  1e-07    2e+00  67 5.0 20
HD103095         -9e+34 -30 5.0 20
''' # noqa: W291
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
    lines = lines[i_secs[-6]:]  # Select Byte-By-Byte section and later lines.
    assert lines == exp_output.splitlines()


exp_coord_cols_output = dict(generic = '''\
--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  8  A8     ---    names   Description of names              
 10- 14  E5.1   ---    e       [-3160000.0/0.01] Description of e
 16- 23  F8.5   ---    d       [22.25/27.25] Description of d    
 25- 31  E7.1   ---    s       [-9e+34/2.0] Description of s     
 33- 35  I3     ---    i       [-30/67] Description of i         
 37- 39  F3.1   ---    sameF   [5.0/5.0] Description of sameF    
 41- 42  I2     ---    sameI   [20] Description of sameI         
 44- 47  F4.1   h      RAh     Right Ascension (hour)            
 49- 52  F4.1   min    RAm     Right Ascension (minute)          
 54- 71  F18.15 s      RAs     Right Ascension (second)          
     73  A1     ---    DE-     Sign of Declination               
 74- 77  F5.1   deg    DEd     Declination (degree)              
 79- 82  F4.1   arcmin DEm     Declination (arcmin)              
 84-101  F18.15 arcsec DEs     Declination (arcsec)              
--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by 1st author ?  / astropy.io.ascii )
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 22.0  2.0 15.450000000007265 -61.0 39.0 34.599996000000601
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 12.0 48.0 15.224407200004890  17.0 46.0 26.496624000004374
''',  # noqa: W291

positive_de = '''\
--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  8  A8     ---    names   Description of names              
 10- 14  E5.1   ---    e       [-3160000.0/0.01] Description of e
 16- 23  F8.5   ---    d       [22.25/27.25] Description of d    
 25- 31  E7.1   ---    s       [-9e+34/2.0] Description of s     
 33- 35  I3     ---    i       [-30/67] Description of i         
 37- 39  F3.1   ---    sameF   [5.0/5.0] Description of sameF    
 41- 42  I2     ---    sameI   [20] Description of sameI         
 44- 47  F4.1   h      RAh     Right Ascension (hour)            
 49- 52  F4.1   min    RAm     Right Ascension (minute)          
 54- 70  F17.14 s      RAs     Right Ascension (second)          
 72- 75  F4.1   deg    DEd     Declination (degree)              
 77- 80  F4.1   arcmin DEm     Declination (arcmin)              
 82- 99  F18.15 arcsec DEs     Declination (arcsec)              
--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by 1st author ?  / astropy.io.ascii )
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 12.0 48.0 15.22440720000489 17.0 46.0 26.496624000004374
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 12.0 48.0 15.22440720000489 17.0 46.0 26.496624000004374
''',  # noqa: W291

galactic = '''\
--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 8  A8     ---    names   Description of names              
10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e
16-23  F8.5   ---    d       [22.25/27.25] Description of d    
25-31  E7.1   ---    s       [-9e+34/2.0] Description of s     
33-35  I3     ---    i       [-30/67] Description of i         
37-39  F3.1   ---    sameF   [5.0/5.0] Description of sameF    
41-42  I2     ---    sameI   [20] Description of sameI         
44-60  F17.13 deg    GLON    Galactic Longitude                
62-79  F18.14 deg    GLAT    Galactic Latitude                 
--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by 1st author ?  / astropy.io.ascii )
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 330.0716395916897 -45.54808048460931
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 330.0716395916897 -45.54808048460931
''',  # noqa: W291

ecliptic = '''\
--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 8  A8     ---    names   Description of names                       
10-14  E5.1   ---    e       [-3160000.0/0.01] Description of e         
16-23  F8.5   ---    d       [22.25/27.25] Description of d             
25-31  E7.1   ---    s       [-9e+34/2.0] Description of s              
33-35  I3     ---    i       [-30/67] Description of i                  
37-39  F3.1   ---    sameF   [5.0/5.0] Description of sameF             
41-42  I2     ---    sameI   [20] Description of sameI                  
44-60  F17.13 deg    ELON    Ecliptic Longitude (geocentrictrueecliptic)
62-79  F18.14 deg    ELAT    Ecliptic Latitude (geocentrictrueecliptic) 
--------------------------------------------------------------------------------

See also:


Acknowledgements:

References:
================================================================================
     (prepared by 1st author ?  / astropy.io.ascii )
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 306.2242086500961 -45.62178985082456
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 306.2242086500961 -45.62178985082456
'''  # noqa: W291
)


def test_write_coord_cols():
    """
    There can only be one such coordinate column in a single table,
    because division of columns into individual component columns requires
    iterating over the table columns, which will have to be done again
    if additional such coordinate columns are present.
    """
    from astropy.table import Column
    from astropy.coordinates import SkyCoord
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name='sameF')
    t.add_column([20, 20], name='sameI')

    # Coordinates of ASASSN-15lh
    coord = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    # Coordinates of ASASSN-14li
    coordp = SkyCoord(192.06343503, 17.77402684, unit=u.deg)
    cols = [Column([coord, coordp]), # Generic coordinate column
            coordp,          # Coordinate column with positive DEC
            coord.galactic,  # Galactic coordinates
            coord.geocentrictrueecliptic  # Ecliptic coordinates
            ]

    # Loop through different types of coordinate columns.
    for col, exp_output in zip(cols, exp_coord_cols_output.values()):
        t['coord'] = col
        out = StringIO()
        t.write(out, format='ascii.cds')
        lines = out.getvalue().splitlines()
        i_secs = [i for i, s in enumerate(lines)
                  if s.startswith(('------', '======='))]
        lines = lines[i_secs[-6]:]  # Select Byte-By-Byte section and later lines.
        # Check the written table.
        assert lines == exp_output.splitlines()

        # Check if the original table columns remains unmodified.
        assert t.colnames == ['names', 'e', 'd', 's', 'i', 'sameF', 'sameI', 'coord']


def test_write_byte_by_byte_bytes_col_format():
    """
    Tests the alignment of Byte counts with respect to hyphen
    in the Bytes column of Byte-By-Byte. The whitespace around the
    hyphen is govered by the number of digits in the total Byte
    count. Single Byte columns should have a single Byte count
    without the hyphen.
    """
    exp_output = '''\
--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  8  A8     ---    names         Description of names              
 10- 14  E5.1   ---    e             [-3160000.0/0.01] Description of e
 16- 23  F8.5   ---    d             [22.25/27.25] Description of d    
 25- 31  E7.1   ---    s             [-9e+34/2.0] Description of s     
 33- 35  I3     ---    i             [-30/67] Description of i         
 37- 39  F3.1   ---    sameF         [5.0/5.0] Description of sameF    
 41- 42  I2     ---    sameI         [20] Description of sameI         
     44  I1     ---    singleByteCol [2] Description of singleByteCol  
 46- 49  F4.1   h      RAh           Right Ascension (hour)            
 51- 53  F3.1   min    RAm           Right Ascension (minute)          
 55- 72  F18.15 s      RAs           Right Ascension (second)          
     74  A1     ---    DE-           Sign of Declination               
 75- 78  F5.1   deg    DEd           Declination (degree)              
 80- 83  F4.1   arcmin DEm           Declination (arcmin)              
 85-100  F16.13 arcsec DEs           Declination (arcsec)              
--------------------------------------------------------------------------------
''' # noqa: W291
    from astropy.coordinates import SkyCoord
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name='sameF')
    t.add_column([20, 20], name='sameI')
    t['coord'] = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    t['singleByteCol'] = [2, 2]
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines)
              if s.startswith(('------', '======='))]
    # Select only the Byte-By-Byte section.
    lines = lines[i_secs[-6]:i_secs[-3]]
    lines.append(lines[0])   # Append a separator line.
    assert lines == exp_output.splitlines()


def test_write_byte_by_byte_wrapping():
    """
    Test line wrapping in the description column of the
    Byte-By-Byte section of the ReadMe.
    """
    exp_output = '''\
--------------------------------------------------------------------------------
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 8  A8     ---    thisIsALongColumnLabel This is a tediously long
                                           description. But they do sometimes
                                           have them. Better to put extra
                                           details in the notes. This is a
                                           tediously long description. But they
                                           do sometimes have them. Better to put
                                           extra details in the notes.
10-14  E5.1   ---    e                      [-3160000.0/0.01] Description of e
16-23  F8.5   ---    d                      [22.25/27.25] Description of d
--------------------------------------------------------------------------------
''' # noqa: W291
    t = ascii.read(test_dat)
    t.remove_columns(['s', 'i'])
    description = 'This is a tediously long description.' \
                  + ' But they do sometimes have them.' \
                  + ' Better to put extra details in the notes. '
    t['names'].description = description * 2
    t['names'].name = 'thisIsALongColumnLabel'
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines)
              if s.startswith(('------', '======='))]
    # Select only the Byte-By-Byte section.
    lines = lines[i_secs[-6]:i_secs[-3]]
    lines.append(lines[0])   # Append a separator line.
    assert lines == exp_output.splitlines()
