# Licensed under a 3-clause BSD style license - see LICENSE.rst
# this test doesn't actually use any online data, it should just be skipped
# by run_tests because it has the remote_data decorator.
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..helper import remote_data
from ..helper import pytest
from ...utils.data import get_pkg_data_filename, download_file


@remote_data
def test_skip_remote_data(pytestconfig):

    # astropy.test() has remote_data=none or remote_data=astropy but we still
    # got here somehow, so fail with a helpful message

    if pytestconfig.getoption('remote_data') == 'none':
        pytest.fail('@remote_data was not skipped with remote_data=none')
    elif pytestconfig.getoption('remote_data') == 'astropy':
        pytest.fail('@remote_data was not skipped with remote_data=astropy')

    # Test Astropy URL
    get_pkg_data_filename('galactic_center/gc_2mass_k.fits')

    # Test non-Astropy URL
    download_file('http://www.google.com')


@remote_data(source='astropy')
def test_skip_remote_data_astropy(pytestconfig):

    # astropy.test() has remote_data=none but we still got here somehow,
    # so fail with a helpful message

    if pytestconfig.getoption('remote_data') == 'none':
        pytest.fail('@remote_data was not skipped with remote_data=none')

    # Test Astropy URL
    get_pkg_data_filename('galactic_center/gc_2mass_k.fits')

    # Test non-Astropy URL
    if pytestconfig.getoption('remote_data') == 'astropy':
        with pytest.raises(Exception) as exc:
            download_file('http://www.google.com')
        assert "An attempt was made to connect to the internet" in str(exc.value)
    else:
        download_file('http://www.google.com')
