# -*- coding: utf-8 -*-

import os

from ..scripts import showtable

ROOT = os.path.abspath(os.path.dirname(__file__))
ASCII_ROOT = os.path.join(ROOT, '..', '..', 'io', 'ascii', 'tests')
FITS_ROOT = os.path.join(ROOT, '..', '..', 'io', 'fits', 'tests')


def test_missing_file(capsys):
    showtable.main(['foobar.fits'])
    out, err = capsys.readouterr()
    assert err.startswith("ERROR: [Errno 2] No such file or directory: "
                          "'foobar.fits'")


def test_fits(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/table.fits')])
    out, err = capsys.readouterr()
    assert out == (' target V_mag\n'
                   '             \n'
                   '------- -----\n'
                   'NGC1001  11.1\n'
                   'NGC1002  12.3\n'
                   'NGC1003  15.2\n')


def test_fits_hdu(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/zerowidth.fits'),
                    '--hdu', 'AIPS OF'])
    out, err = capsys.readouterr()
    assert out.startswith(
        '  TIME   SOURCE ID ANTENNA NO. SUBARRAY FREQ ID ANT FLAG STATUS 1\n'
        '  DAYS                                                           \n'
        '-------- --------- ----------- -------- ------- -------- --------\n'
        '0.144387         1          10        1       1        4        4\n'
    )


def test_csv(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/simple_csv.csv')])
    out, err = capsys.readouterr()
    assert out == (' a   b   c \n'
                   '           \n'
                   '--- --- ---\n'
                   '  1   2   3\n'
                   '  4   5   6\n')


def test_ascii_format(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/commented_header.dat'),
                    '--format', 'ascii.commented_header'])
    out, err = capsys.readouterr()
    assert out == (' a   b   c \n'
                   '           \n'
                   '--- --- ---\n'
                   '  1   2   3\n'
                   '  4   5   6\n')
