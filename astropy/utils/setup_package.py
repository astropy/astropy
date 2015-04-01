# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import

from distutils.core import Extension
from os.path import dirname, join, relpath

ASTROPY_UTILS_ROOT = dirname(__file__)


def get_extensions():
    return [
        Extension('astropy.utils._compiler',
                  [relpath(join(ASTROPY_UTILS_ROOT, 'src', 'compiler.c'))])
    ]


def get_package_data():
    # Installs the testing data files
    return {
        'astropy.utils': ['data/*.txt'],
        'astropy.utils.tests': [
            'data/test_package/*.py',
            'data/test_package/data/*.txt',
            'data/*.dat',
            'data/*.dat.gz',
            'data/*.dat.bz2',
            'data/*.txt',
            'data/.hidden_file.txt',
            'data/*.cfg'],
        'astropy.utils.iers': [
            'data/ReadMe.eopc04_IAU2000',
            'data/ReadMe.finals2000A',
            'data/eopc04_IAU2000.62-now']
    }


def requires_2to3():
    return False
