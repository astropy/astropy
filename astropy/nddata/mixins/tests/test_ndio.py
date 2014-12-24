from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..nddata import NDData
from ..io import NDIOMixin


# Define minimal class that uses the I/O mixin
class NDDataIO(NDIOMixin, NDData):
    pass


def test_simple_write_read(tmpdir):
    ndd = NDDataIO([1, 2, 3])
    assert hasattr(ndd, 'read')
    assert hasattr(ndd, 'write')
