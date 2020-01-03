# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import pytest
import numpy as np

from astropy.io.ascii import read

ROOT = os.path.abspath(os.path.dirname(__file__))


try:
    import bz2  # noqa
except ImportError:
    HAS_BZ2 = False
else:
    HAS_BZ2 = True

try:
    import lzma  # noqa
except ImportError:
    HAS_XZ = False
else:
    HAS_XZ = True


@pytest.mark.parametrize('filename', ['data/daophot.dat.gz', 'data/latex1.tex.gz',
                                      'data/short.rdb.gz'])
def test_gzip(filename):
    t_comp = read(os.path.join(ROOT, filename))
    t_uncomp = read(os.path.join(ROOT, filename.replace('.gz', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())


@pytest.mark.xfail('not HAS_BZ2')
@pytest.mark.parametrize('filename', ['data/short.rdb.bz2', 'data/ipac.dat.bz2'])
def test_bzip2(filename):
    t_comp = read(os.path.join(ROOT, filename))
    t_uncomp = read(os.path.join(ROOT, filename.replace('.bz2', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())


@pytest.mark.xfail('not HAS_XZ')
@pytest.mark.parametrize('filename', ['data/short.rdb.xz', 'data/ipac.dat.xz'])
def test_xz(filename):
    t_comp = read(os.path.join(ROOT, filename))
    t_uncomp = read(os.path.join(ROOT, filename.replace('.xz', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())
