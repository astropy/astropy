# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest
import numpy as np

from astropy.io.ascii import read
from astropy.utils.data import get_pkg_data_filename

# NOTE: Python can be built without bz2 or lzma
from astropy.utils.compat.optional_deps import HAS_BZ2, HAS_LZMA  # noqa


@pytest.mark.parametrize('filename', ['data/daophot.dat.gz', 'data/latex1.tex.gz',
                                      'data/short.rdb.gz'])
def test_gzip(filename):
    t_comp = read(get_pkg_data_filename(filename))
    t_uncomp = read(get_pkg_data_filename(filename.replace('.gz', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())


@pytest.mark.xfail('not HAS_BZ2')
@pytest.mark.parametrize('filename', ['data/short.rdb.bz2', 'data/ipac.dat.bz2'])
def test_bzip2(filename):
    t_comp = read(get_pkg_data_filename(filename))
    t_uncomp = read(get_pkg_data_filename(filename.replace('.bz2', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())


@pytest.mark.xfail('not HAS_LZMA')
@pytest.mark.parametrize('filename', ['data/short.rdb.xz', 'data/ipac.dat.xz'])
def test_xz(filename):
    t_comp = read(get_pkg_data_filename(filename))
    t_uncomp = read(get_pkg_data_filename(filename.replace('.xz', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())
