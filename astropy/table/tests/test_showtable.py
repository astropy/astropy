# -*- coding: utf-8 -*-

import os

from ..scripts import showtable

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
    assert out == ('<Table length=3>{0}'
                   ' name   dtype {0}'
                   '------ -------{0}'
                   'target   str20{0}'
                   ' V_mag float32{0}'
                   '\n').format(os.linesep)


def test_stats(capsys):
    showtable.main([os.path.join(FITS_ROOT, 'data/table.fits'), '--stats'])
    out, err = capsys.readouterr()
    assert out == ('<Table length=3>{0}'
                   ' name    mean    std   min  max {0}'
                   '------ ------- ------- ---- ----{0}'
                   'target      --      --   --   --{0}'
                   ' V_mag 12.8667 1.72111 11.1 15.2{0}').format(os.linesep)


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


def test_ascii_delimiter(capsys):
    showtable.main([os.path.join(ASCII_ROOT, 't/simple2.txt'),
                    '--format', 'ascii', '--delimiter', '|'])
    out, err = capsys.readouterr()
    assert out == (
        "obsid redshift  X    Y      object   rad \n"
        "                                         \n"
        "----- -------- ---- ---- ----------- ----\n"
        " 3102     0.32 4167 4085 Q1250+568-A  9.0\n"
        " 3102     0.32 4706 3916 Q1250+568-B 14.0\n"
        "  877     0.22 4378 3892 'Source 82' 12.5\n"
    )


def test_votable(capsys):
    showtable.main([os.path.join(VOTABLE_ROOT, 'data/regression.xml'),
                    '--table-id', 'main_table', '--max-width', '50'])
    out, err = capsys.readouterr()
    assert out == (
        '   string_test    string_test_2 ... bitarray2 [16]\n'
        '                                ...               \n'
        '----------------- ------------- ... --------------\n'
        '    String & test    Fixed stri ...  True .. False\n'
        'String &amp; test    0123456789 ...       -- .. --\n'
        '             XXXX          XXXX ...       -- .. --\n'
        '                                ...       -- .. --\n'
        '                                ...       -- .. --\n'
    )
