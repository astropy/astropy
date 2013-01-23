# Licensed under a 3-clause BSD style license - see LICENSE.rst
import io
import os

from ...tests.helper import pytest
from ..compat import gzip

pytestmark = pytest.mark.skipif("sys.version_info < (3,0)")


def test_gzip(tmpdir):
    fd = gzip.GzipFile(str(tmpdir.join("test.gz")), 'wb')
    fd = io.TextIOWrapper(fd, encoding='utf8')


def test_gzip2(tmpdir):
    with gzip.GzipFile(str(tmpdir.join("test.gz")), 'wb') as fd:
        pass

