# Licensed under a 3-clause BSD style license - see LICENSE.rst


from os.path import join

import numpy as np
from setuptools import Extension


def get_extensions(build_type="release"):
    VO_DIR = "astropy/io/votable/src"

    return [
        Extension(
            "astropy.io.votable.tablewriter",
            [join(VO_DIR, "tablewriter.c")],
            include_dirs=[VO_DIR],
        ),
        Extension(
            "astropy.io.votable.fast_converters",
            [join(VO_DIR, "fast_converters.pyx")],
            include_dirs=[np.get_include(), VO_DIR],
        ),
    ]
