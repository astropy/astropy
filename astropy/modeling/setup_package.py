# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {
        'astropy.modeling.tests': ['data/*.fits', 'data/*.hdr',
                                   '../../wcs/tests/maps/*.hdr']
    }


def requires_2to3():
    return False
