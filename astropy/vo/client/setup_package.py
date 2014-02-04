# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {
        'astropy.vo.client.tests': ['data/*.xml', 'data/*.json']}


def requires_2to3():
    return False
