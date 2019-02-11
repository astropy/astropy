# Licensed under a 3-clause BSD style license

import os
from distutils.extension import Extension

ROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    sources = [os.path.join(ROOT, 'cparser.pyx'),
               os.path.join(ROOT, 'src', 'tokenizer.c')]
    ascii_ext = Extension(
        name="astropy.io.ascii.cparser",
        include_dirs=["numpy"],
        sources=sources)
    return [ascii_ext]


def get_package_data():
    # Installs the testing data files.  Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {
        'astropy.io.ascii.tests': ['data/vizier/ReadMe',
                                   'data/vizier/table1.dat',
                                   'data/vizier/table5.dat',
                                   'data/apostrophe.rdb',
                                   'data/apostrophe.tab',
                                   'data/bad.txt',
                                   'data/bars_at_ends.txt',
                                   'data/cds.dat',
                                   'data/cds_malformed.dat',
                                   'data/cds/glob/ReadMe',
                                   'data/cds/glob/lmxbrefs.dat',
                                   'data/cds/multi/ReadMe',
                                   'data/cds/multi/lhs2065.dat',
                                   'data/cds/multi/lp944-20.dat',
                                   'data/cds2.dat',
                                   'data/commented_header.dat',
                                   'data/commented_header2.dat',
                                   'data/continuation.dat',
                                   'data/daophot.dat',
                                   'data/daophot2.dat',
                                   'data/daophot3.dat',
                                   'data/daophot4.dat',
                                   'data/sextractor.dat',
                                   'data/sextractor2.dat',
                                   'data/sextractor3.dat',
                                   'data/daophot.dat.gz',
                                   'data/fill_values.txt',
                                   'data/html.html',
                                   'data/html2.html',
                                   'data/ipac.dat',
                                   'data/ipac.dat.bz2',
                                   'data/ipac.dat.xz',
                                   'data/latex1.tex',
                                   'data/latex1.tex.gz',
                                   'data/latex2.tex',
                                   'data/latex3.tex',
                                   'data/nls1_stackinfo.dbout',
                                   'data/no_data_cds.dat',
                                   'data/no_data_daophot.dat',
                                   'data/no_data_sextractor.dat',
                                   'data/no_data_ipac.dat',
                                   'data/no_data_with_header.dat',
                                   'data/no_data_without_header.dat',
                                   'data/short.rdb',
                                   'data/short.rdb.bz2',
                                   'data/short.rdb.gz',
                                   'data/short.rdb.xz',
                                   'data/short.tab',
                                   'data/simple.txt',
                                   'data/simple2.txt',
                                   'data/simple3.txt',
                                   'data/simple4.txt',
                                   'data/simple5.txt',
                                   'data/space_delim_blank_lines.txt',
                                   'data/space_delim_no_header.dat',
                                   'data/space_delim_no_names.dat',
                                   'data/test4.dat',
                                   'data/test5.dat',
                                   'data/vots_spec.dat',
                                   'data/whitespace.dat',
                                   'data/simple_csv.csv',
                                   'data/simple_csv_missing.csv',
                                   'data/fixed_width_2_line.txt',
                                   'data/cds/description/ReadMe',
                                   'data/cds/description/table.dat',
                                   ]
    }
