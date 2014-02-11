# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {
        str('astropy.config.tests'): ['data/*.cfg']
    }


def requires_2to3():
    return False
