# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {'astropy.coordinates.tests.accuracy': ['*.csv'],
            'astropy.coordinates':['data/*.dat']}


def requires_2to3():
    return False
