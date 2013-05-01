# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {
        'astropy.vo.validator': ['data/*.txt'],
        'astropy.vo.validator.tests': ['data/*.json', 'data/*.xml',
                                       'data/*.out']}
