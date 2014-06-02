# **Please Note**: ``astropy.sphinx`` exists only for backward-compatibility
# purposes - it has now been moved to the separate astropy-helpers package,
# located at https://github.com/astropy/astropy-helpers. Any new development or
# bug fixes should be done there.
# Licensed under a 3-clause BSD style license - see LICENSE.rst

def get_package_data():
    # Install the theme files
    return {
        'astropy.sphinx': [
            'ext/templates/*/*',
            'themes/bootstrap-astropy/*.*',
            'themes/bootstrap-astropy/static/*.*']}
