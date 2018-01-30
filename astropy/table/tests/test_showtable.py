import os
import re

from ..scripts import showtable
from ...utils.compat import NUMPY_LT_1_14

ROOT = os.path.abspath(os.path.dirname(__file__))
ASCII_ROOT = os.path.join(ROOT, '..', '..', 'io', 'ascii', 'tests')
FITS_ROOT = os.path.join(ROOT, '..', '..', 'io', 'fits', 'tests')
VOTABLE_ROOT = os.path.join(ROOT, '..', '..', 'io', 'votable', 'tests')


def test_missing_file(capsys):
    showtable.main(['foobar.fits'])
    out, err = capsys.readouterr()
    assert err.startswith("ERROR: [Errno 2] No such file or directory: "
                          "'foobar.fits'")


def test_info(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/table.fits'), '--info'])
    out, err = capsys.readouterr()
    assert out.splitlines() == ['<Table length=3>',
                                ' name   dtype ',
                                '------ -------',
                                'target bytes20',
                                ' V_mag float32']


def test_stats(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/table.fits'), '--stats'])
    out, err = capsys.readouterr()
    if NUMPY_LT_1_14:
        expected = ['<Table length=3>',
                    ' name    mean    std   min  max ',
                    '------ ------- ------- ---- ----',
                    'target      --      --   --   --',
                    ' V_mag 12.8667 1.72111 11.1 15.2']
    else:
        expected = ['<Table length=3>',
                    ' name     mean      std    min  max ',
                    '------ --------- --------- ---- ----',
                    'target        --        --   --   --',
                    ' V_mag 12.86666[0-9]? 1.7211105 11.1 15.2']

    out = out.splitlines()
    assert out[:4] == expected[:4]
    # Here we use re.match as in some cases one of the values above is
    # platform-dependent.
    assert re.match(expected[4], out[4]) is not None


def test_fits(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/table.fits')])
    out, err = capsys.readouterr()
    assert out.splitlines() == [' target V_mag',
                                '------- -----',
                                'NGC1001  11.1',
                                'NGC1002  12.3',
                                'NGC1003  15.2']


def test_fits_hdu(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/zerowidth.fits'),
                    '--hdu', 'AIPS OF'])
    out, err = capsys.readouterr()
    if NUMPY_LT_1_14:
        assert out.startswith(
            '  TIME   SOURCE ID ANTENNA NO. SUBARRAY FREQ ID ANT FLAG STATUS 1\n'
            '  DAYS                                                           \n'
            '-------- --------- ----------- -------- ------- -------- --------\n'
            '0.144387         1          10        1       1        4        4\n')
    else:
        assert out.startswith(
            '   TIME    SOURCE ID ANTENNA NO. SUBARRAY FREQ ID ANT FLAG STATUS 1\n'
            '   DAYS                                                            \n'
            '---------- --------- ----------- -------- ------- -------- --------\n'
            '0.14438657         1          10        1       1        4        4\n')


def test_csv(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/simple_csv.csv')])
    out, err = capsys.readouterr()
    assert out.splitlines() == [' a   b   c ',
                                '--- --- ---',
                                '  1   2   3',
                                '  4   5   6']


def test_ascii_format(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/commented_header.dat'),
                    '--format', 'ascii.commented_header'])
    out, err = capsys.readouterr()
    assert out.splitlines() == [' a   b   c ',
                                '--- --- ---',
                                '  1   2   3',
                                '  4   5   6']


def test_ascii_delimiter(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/simple2.txt'),
                    '--format', 'ascii', '--delimiter', '|'])
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        "obsid redshift  X    Y      object   rad ",
        "----- -------- ---- ---- ----------- ----",
        " 3102     0.32 4167 4085 Q1250+568-A  9.0",
        " 3102     0.32 4706 3916 Q1250+568-B 14.0",
        "  877     0.22 4378 3892 'Source 82' 12.5",
    ]


def test_votable(capsys):
    showtable.main([os.path.join(VOTABLE_ROOT, 'data/regression.xml'),
                    '--table-id', 'main_table', '--max-width', '50'])
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        '   string_test    string_test_2 ... bitarray2 [16]',
        '----------------- ------------- ... --------------',
        '    String & test    Fixed stri ...  True .. False',
        'String &amp; test    0123456789 ...       -- .. --',
        '             XXXX          XXXX ...       -- .. --',
        '                                ...       -- .. --',
        '                                ...       -- .. --',
    ]


def test_max_lines(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/cds2.dat'),
                    '--format', 'ascii.cds', '--max-lines', '7',
                    '--max-width', '30'])
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        '      SST       ... Note',
        '                ...     ',
        '--------------- ... ----',
        '041314.1+281910 ...   --',
        '            ... ...  ...',
        '044427.1+251216 ...   --',
        '044642.6+245903 ...   --',
        'Length = 215 rows',
    ]


def test_show_dtype(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/table.fits'),
                    '--show-dtype'])
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        ' target  V_mag ',
        'bytes20 float32',
        '------- -------',
        'NGC1001    11.1',
        'NGC1002    12.3',
        'NGC1003    15.2',
    ]


def test_hide_unit(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/cds.dat'),
                    '--format', 'ascii.cds'])
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        'Index RAh RAm  RAs  DE- DEd  DEm    DEs   Match Class  AK Fit ',
        '       h  min   s       deg arcmin arcsec             mag     ',
        '----- --- --- ----- --- --- ------ ------ ----- ----- --- ----',
        '    1   3  28 39.09   +  31      6    1.9    --    I*  -- 1.35',
    ]

    showtable.main([os.path.join(ASCII_ROOT, 't/cds.dat'),
                    '--format', 'ascii.cds', '--hide-unit'])
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        'Index RAh RAm  RAs  DE- DEd DEm DEs Match Class  AK Fit ',
        '----- --- --- ----- --- --- --- --- ----- ----- --- ----',
        '    1   3  28 39.09   +  31   6 1.9    --    I*  -- 1.35',
    ]
