# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os


def get_package_data():
    return {
            'astropy.samp': [os.path.join('data', '*')],
            'astropy.samp.tests': [os.path.join('data', '*')]
           }
