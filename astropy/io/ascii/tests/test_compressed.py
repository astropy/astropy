# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import sys
import numpy as np

from ....tests.helper import pytest
from .. import read

ROOT = os.path.abspath(os.path.dirname(__file__))


try:
    import bz2
except ImportError:
    HAS_BZ2 = False
else:
    HAS_BZ2 = True

try:
    if sys.version_info >= (3,3,0):
        import lzma
    else:
        from backports import lzma
except ImportError:
    HAS_XZ = False
else:
    HAS_XZ = True


@pytest.mark.parametrize('filename', ['t/daophot.dat.gz', 't/latex1.tex.gz',
                                      't/short.rdb.gz'])
def test_gzip(filename):
    t_comp = read(os.path.join(ROOT, filename))
    t_uncomp = read(os.path.join(ROOT, filename.replace('.gz', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())


@pytest.mark.xfail('not HAS_BZ2')
@pytest.mark.parametrize('filename', ['t/short.rdb.bz2', 't/ipac.dat.bz2'])
def test_bzip2(filename):
    t_comp = read(os.path.join(ROOT, filename))
    t_uncomp = read(os.path.join(ROOT, filename.replace('.bz2', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())


@pytest.mark.xfail('not HAS_XZ')
@pytest.mark.parametrize('filename', ['t/short.rdb.xz', 't/ipac.dat.xz'])
def test_xz(filename):
    t_comp = read(os.path.join(ROOT, filename))
    t_uncomp = read(os.path.join(ROOT, filename.replace('.xz', '')))
    assert t_comp.dtype.names == t_uncomp.dtype.names
    assert np.all(t_comp.as_array() == t_uncomp.as_array())
