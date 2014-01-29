# Licensed under a 3-clause BSD style license - see LICENSE.rst


import os


def get_package_data():
    paths = [os.path.join('js', '*.js')]
    return {'astropy.extern': paths}


def requires_2to3():
    return False
