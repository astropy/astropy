# Licensed under a 3-clause BSD style license


def get_package_data():
    # Installs the testing data files.  Unable to get package_data
    # to deal with a directory hierarchy of files, so just explicitly list.
    return {
        'astropy.io.ascii.tests': ['t/vizier/ReadMe',
                                   't/vizier/table1.dat',
                                   't/vizier/table5.dat',
                                   't/apostrophe.rdb',
                                   't/apostrophe.tab',
                                   't/bad.txt',
                                   't/bars_at_ends.txt',
                                   't/cds.dat',
                                   't/cds/glob/ReadMe',
                                   't/cds/glob/lmxbrefs.dat',
                                   't/cds/multi/ReadMe',
                                   't/cds/multi/lhs2065.dat',
                                   't/cds2.dat',
                                   't/commented_header.dat',
                                   't/continuation.dat',
                                   't/daophot.dat',
                                   't/fill_values.txt',
                                   't/ipac.dat',
                                   't/latex1.tex',
                                   't/latex2.tex',
                                   't/nls1_stackinfo.dbout',
                                   't/no_data_cds.dat',
                                   't/no_data_daophot.dat',
                                   't/no_data_ipac.dat',
                                   't/no_data_with_header.dat',
                                   't/no_data_without_header.dat',
                                   't/short.rdb',
                                   't/short.tab',
                                   't/simple.txt',
                                   't/simple2.txt',
                                   't/simple3.txt',
                                   't/simple4.txt',
                                   't/simple5.txt',
                                   't/space_delim_blank_lines.txt',
                                   't/space_delim_no_header.dat',
                                   't/space_delim_no_names.dat',
                                   't/test4.dat',
                                   't/test5.dat',
                                   't/vots_spec.dat',
                                   't/whitespace.dat',
                                   ]
        }
