import io
import os

from ...tests.helper import pytest
from .. import gzip

def test_gzip(tmpdir):
    fd = gzip.GzipFile(str(tmpdir.join("test.gz")), 'wb')
    fd = io.TextIOWrapper(fd, encoding='utf8')
