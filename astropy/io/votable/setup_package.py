# Licensed under a 3-clause BSD style license - see LICENSE.rst

from distutils.core import Extension
from os.path import join


def get_extensions(build_type='release'):
    VO_DIR = 'astropy/io/votable/src'

    return [Extension(
        "astropy.io.votable.tablewriter",
        [join(VO_DIR, "tablewriter.c")],
        include_dirs=[VO_DIR])]
