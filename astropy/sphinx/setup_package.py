# Licensed under a 3-clause BSD style license - see LICENSE.rst

def get_package_data():
    # Install the theme files
    return {
        'astropy.sphinx': [
            'ext/templates/*/*',
            'themes/bootstrap-astropy/*.*',
            'themes/bootstrap-astropy/static/*.*']}
