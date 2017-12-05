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
    assert out.splitlines() == ['<Table length=3>',
                                ' name   dtype ',
                                '------ -------',
                                'target   str20',
                                ' V_mag float32']


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
    assert out.splitlines() == [' target V_mag',
                                '------- -----',
                                'NGC1001  11.1',
                                'NGC1002  12.3',
                                'NGC1003  15.2']


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
