# Licensed under a 3-clause BSD style license - see LICENSE.rst

from os.path import join

from setuptools import Extension


def get_extensions(build_type="release"):
    VO_DIR = "astropy/io/votable/src"

    return [
        Extension(
            "astropy.io.votable.tablewriter",
            [join(VO_DIR, "tablewriter.c")],
            include_dirs=[VO_DIR],
        )
    ]
