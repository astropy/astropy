# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {
        'astropy.tests': ['coveragerc'],
        'astropy.tests.tests': ['data/open_file_detection.txt']}


def requires_2to3():
    return False
