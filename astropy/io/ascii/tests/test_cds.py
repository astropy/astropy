# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This module tests some methods related to ``CDS`` format
reader/writer.
Requires `pyyaml <https://pyyaml.org/>`_ to be installed.
"""
import os
from io import StringIO

from astropy.io import ascii
from astropy import units as u
from astropy.table import Table
from astropy.table import Column, MaskedColumn
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename


test_dat = ['names e d s i',
            'HD81809 1E-7 22.25608 +2 67',
            'HD103095 -31.6e5 +27.2500 -9E34 -30']


def test_roundtrip_mrt_table():
    """
    Tests whether or not the CDS writer can roundtrip a table,
    i.e. read a table to ``Table`` object and write it exactly
    as it is back to a file. Since, presently CDS uses a
    MRT format template while writing, only the Byte-By-Byte
    and the data section of the table can be compared between
    original and the newly written table.

    Further, the CDS Reader does not have capability to recognize
    column format from the header of a CDS/MRT table, so this test
    can work for a limited set of simple tables, which don't have
    whitespaces in the column values or mix-in columns. Because of
    this the written table output cannot be directly matched with
    the original file and have to be checked against a list of lines.
    Masked columns are read properly though, and thus are being tested
    during round-tripping.

    The difference between ``cdsFunctional2.dat`` file and ``exp_output``
    is the following:
        * Metadata is different because MRT template is used for writing.
        * Spacing between ``Label`` and ``Explanations`` column in the
            Byte-By-Byte.
        * Units are written as ``[cm.s-2]`` and not ``[cm/s2]``, since both
            are valid according to CDS/MRT standard.
    """
    exp_output = '''\
================================================================================
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 7  A7       ---    ID       Star ID                              
 9-12  I4       K      Teff     [4337/4654] Effective temperature    
14-17  F4.2   [cm.s-2] logg     [0.77/1.28] Surface gravity          
19-22  F4.2     km.s-1 vturb    [1.23/1.82] Micro-turbulence velocity
24-28  F5.2     [-]    [Fe/H]   [-2.11/-1.5] Metallicity             
30-33  F4.2     [-]    e_[Fe/H] ? rms uncertainty on [Fe/H]          
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
S05-5   4337 0.77 1.80 -2.07     
S08-229 4625 1.23 1.23 -1.50     
S05-10  4342 0.91 1.82 -2.11 0.14
S05-47  4654 1.28 1.74 -1.64 0.16
''' # noqa: W291
    dat = get_pkg_data_filename('data/cdsFunctional2.dat',
                                package='astropy.io.ascii.tests')
    t = Table.read(dat, format='ascii.cds')
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_bbb = lines.index('=' * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
    assert lines == exp_output.splitlines()


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
Title:
Authors:
Table:
================================================================================
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
Notes:
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67
HD103095 -3e+06 27.25000  -9e+34 -30
''' # noqa: W291
    t = ascii.read(test_dat)
    out = StringIO()
    t.write(out, format='ascii.cds')
    assert out.getvalue().splitlines() == exp_output.splitlines()


def test_write_empty_table():
    out = StringIO()
    import pytest
    with pytest.raises(NotImplementedError):
        Table().write(out, format='ascii.cds')


def test_write_null_data_values():
    exp_output = ['HD81809  1e-07  22.25608  2.0e+00  67',
                  'HD103095 -3e+06 27.25000 -9.0e+34 -30',
                  'Sun                       5.3e+27    ']
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
================================================================================
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 8  A8     ---    names   Description of names          
10-14  E5.1   ---    e       [0.0/0.01]? Description of e  
16-17  F2.0   ---    d       ? Description of d            
19-25  E7.1   ---    s       [-9e+34/2.0] Description of s 
27-29  I3     ---    i       [-30/67] Description of i     
31-33  F3.1   ---    sameF   [5.0/5.0] Description of sameF
35-36  I2     ---    sameI   [20] Description of sameI     
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
HD81809  1e-07    2e+00  67 5.0 20
HD103095         -9e+34 -30 5.0 20
''' # noqa: W291
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name='sameF')
    t.add_column([20, 20], name='sameI')
    t['e'] = MaskedColumn(t['e'], mask=[False, True])
    t['d'] = MaskedColumn(t['d'], mask=[True, True])
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_bbb = lines.index('=' * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
    assert lines == exp_output.splitlines()


exp_coord_cols_output = dict(generic = '''\
================================================================================
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
44-47  F4.1   h      RAh     Right Ascension (hour)            
49-52  F4.1   min    RAm     Right Ascension (minute)          
54-68  F15.12 s      RAs     Right Ascension (second)          
   70  A1     ---    DE-     Sign of Declination               
71-74  F5.1   deg    DEd     Declination (degree)              
76-79  F4.1   arcmin DEm     Declination (arcmin)              
81-95  F15.12 arcsec DEs     Declination (arcsec)              
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 22.0  2.0 15.450000000007 -61.0 39.0 34.599996000001
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 12.0 48.0 15.224407200005  17.0 46.0 26.496624000004
''',  # noqa: W291

positive_de = '''\
================================================================================
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
44-47  F4.1   h      RAh     Right Ascension (hour)            
49-52  F4.1   min    RAm     Right Ascension (minute)          
54-68  F15.12 s      RAs     Right Ascension (second)          
70-73  F4.1   deg    DEd     Declination (degree)              
75-78  F4.1   arcmin DEm     Declination (arcmin)              
80-94  F15.12 arcsec DEs     Declination (arcsec)              
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 12.0 48.0 15.224407200005 17.0 46.0 26.496624000004
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 12.0 48.0 15.224407200005 17.0 46.0 26.496624000004
''',  # noqa: W291

galactic = '''\
================================================================================
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
44-59  F16.12 deg    GLON    Galactic Longitude                
61-76  F16.12 deg    GLAT    Galactic Latitude                 
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 330.071639591690 -45.548080484609
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 330.071639591690 -45.548080484609
''',  # noqa: W291

ecliptic = '''\
================================================================================
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
44-59  F16.12 deg    ELON    Ecliptic Longitude (geocentrictrueecliptic)
61-76  F16.12 deg    ELAT    Ecliptic Latitude (geocentrictrueecliptic) 
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
HD81809  1e-07  22.25608   2e+00  67 5.0 20 306.224208650096 -45.621789850825
HD103095 -3e+06 27.25000  -9e+34 -30 5.0 20 306.224208650096 -45.621789850825
'''  # noqa: W291
)


def test_write_coord_cols():
    """
    There can only be one such coordinate column in a single table,
    because division of columns into individual component columns requires
    iterating over the table columns, which will have to be done again
    if additional such coordinate columns are present.
    """
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
        i_bbb = lines.index('=' * 80)
        lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
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
================================================================================
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  8  A8     ---    names         Description of names              
 10- 21  E12.6  ---    e             [-3160000.0/0.01] Description of e
 23- 30  F8.5   ---    d             [22.25/27.25] Description of d    
 32- 38  E7.1   ---    s             [-9e+34/2.0] Description of s     
 40- 42  I3     ---    i             [-30/67] Description of i         
 44- 46  F3.1   ---    sameF         [5.0/5.0] Description of sameF    
 48- 49  I2     ---    sameI         [20] Description of sameI         
     51  I1     ---    singleByteCol [2] Description of singleByteCol  
 53- 56  F4.1   h      RAh           Right Ascension (hour)            
 58- 60  F3.1   min    RAm           Right Ascension (minute)          
 62- 76  F15.12 s      RAs           Right Ascension (second)          
     78  A1     ---    DE-           Sign of Declination               
 79- 82  F5.1   deg    DEd           Declination (degree)              
 84- 87  F4.1   arcmin DEm           Declination (arcmin)              
 89-103  F15.12 arcsec DEs           Declination (arcsec)              
--------------------------------------------------------------------------------
''' # noqa: W291
    t = ascii.read(test_dat)
    t.add_column([5.0, 5.0], name='sameF')
    t.add_column([20, 20], name='sameI')
    t['coord'] = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    t['singleByteCol'] = [2, 2]
    t['e'].format = '.5E'
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_secs = [i for i, s in enumerate(lines)
              if s.startswith(('------', '======='))]
    # Select only the Byte-By-Byte section.
    lines = lines[i_secs[0]:i_secs[-2]]
    lines.append('-'*80)   # Append a separator line.
    assert lines == exp_output.splitlines()


def test_write_byte_by_byte_wrapping():
    """
    Test line wrapping in the description column of the
    Byte-By-Byte section of the ReadMe.
    """
    exp_output = '''\
================================================================================
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
    lines = lines[i_secs[0]:i_secs[-2]]
    lines.append('-'*80)   # Append a separator line.
    print(lines)
    for l, ll in zip(lines, exp_output.splitlines()):
        print(l)
        print(ll)
    assert lines == exp_output.splitlines()


def test_write_mixin_and_broken_cols():
    """
    Tests convertion to string values for ``mix-in`` columns other than
    ``SkyCoord`` and for columns with only partial ``SkyCoord`` values.
    """
    exp_output = '''\
================================================================================
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
  1-  7  A7     ---    name    Description of name   
  9- 74  A66    ---    Unknown Description of Unknown
 76-114  A39    ---    Unknown Description of Unknown
116-138  A23    ---    Unknown Description of Unknown
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
HD81809 <SkyCoord (ICRS): (ra, dec) in deg
    (330.564375, -61.65961111)> (0.41342785, -0.23329341, -0.88014294)  2019-01-01 00:00:00.000
random  12                                                                 (0.41342785, -0.23329341, -0.88014294)  2019-01-01 00:00:00.000
''' # noqa: W291
    t = Table()
    t['name'] = ['HD81809']
    coord = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    t['coord'] = Column(coord)
    t.add_row(['random', 12])
    t['cart'] = coord.cartesian
    t['time'] = Time('2019-1-1')
    out = StringIO()
    t.write(out, format='ascii.cds')
    lines = out.getvalue().splitlines()
    i_bbb = lines.index('=' * 80)
    lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
    # Check the written table.
    assert lines == exp_output.splitlines()


def test_write_extra_SkyCoord_cols():
    """
    Tests output for cases when table contains multiple ``SkyCoord``
    columns.
    """
    exp_output = '''\
================================================================================
Byte-by-byte Description of file: table.dat
--------------------------------------------------------------------------------
 Bytes Format Units  Label     Explanations
--------------------------------------------------------------------------------
 1- 7  A7     ---    name    Description of name     
 9-24  A16    ---    Unknown Description of Unknown  
26-29  F4.1   h      RAh     Right Ascension (hour)  
31-33  F3.1   min    RAm     Right Ascension (minute)
35-52  F18.15 s      RAs     Right Ascension (second)
   54  A1     ---    DE-     Sign of Declination     
55-58  F5.1   deg    DEd     Declination (degree)    
60-63  F4.1   arcmin DEm     Declination (arcmin)    
65-80  F16.13 arcsec DEs     Declination (arcsec)    
--------------------------------------------------------------------------------
Notes:
--------------------------------------------------------------------------------
HD81809 330.564 -61.6596 22.0 2.0 15.450000000007265 -61.0 39.0 34.5999960000006
''' # noqa: W291
    import pytest
    t = Table()
    t['name'] = ['HD81809']
    coord = SkyCoord(330.564375, -61.65961111, unit=u.deg)
    t['coord1'] = coord
    t['coord2'] = coord
    out = StringIO()
    with pytest.raises(UserWarning):
        t.write(out, format='ascii.cds')
        lines = out.getvalue().splitlines()
        i_bbb = lines.index('=' * 80)
        lines = lines[i_bbb:]  # Select Byte-By-Byte section and later lines.
        # Check the written table.
        assert lines == exp_output.splitlines()
