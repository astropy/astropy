# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {'astropy.nddata.tests': ['data/*.fits']}


def requires_2to3():
    return False
