# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import

import os


def get_package_data():
    return {
            'astropy.samp': [os.path.join('data', '*')],
            'astropy.samp.tests': [os.path.join('data', '*')]
           }


def requires_2to3():
    return False
