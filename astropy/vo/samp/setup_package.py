# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os


def get_package_data():
    return {
            'astropy.vo.samp': [os.path.join('data', '*')],
            'astropy.vo.samp.tests': [os.path.join('data', '*')]
           }


def requires_2to3():
    return False
